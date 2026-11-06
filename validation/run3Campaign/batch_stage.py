#!/usr/bin/env python3
"""Execute one frozen Run-3 campaign DAG node in local Condor scratch."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import urlsplit

REPO = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def run(*args):
    subprocess.run(list(map(str, args)), check=True, cwd=REPO)


def load(path):
    return json.loads(Path(path).read_text())


def save(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def xrd_path(path):
    return 'root://eoshome-j.cern.ch/'+str(path).replace('/eos/home-j/', '/eos/user/j/')


def get(url, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    run('xrdcp', '--silent', '--retry', '3', '--cksum', 'adler32', url, path)
    return path


def mkdir_url(url):
    parsed = urlsplit(url)
    run('xrdfs', parsed.scheme+'://'+parsed.netloc, 'mkdir', '-p', '/' + parsed.path.strip('/'))


def publish(directory, url):
    """Publish checked files before completion marker; never infer completion from filenames."""
    directory = Path(directory)
    paths = sorted(p for p in directory.rglob('*') if p.is_file() and p.name != 'stage_complete.json')
    marker = {'status': 'COMPLETE', 'files': {str(p.relative_to(directory)): {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in paths}}
    save(directory/'stage_complete.json', marker)
    paths.append(directory/'stage_complete.json')
    for parent in sorted({str(p.relative_to(directory).parent) for p in paths}):
        mkdir_url(url if parent == '.' else url+'/'+parent)
    for path in paths:
        run('xrdcp', '--silent', '--retry', '3', '--cksum', 'adler32', '--posc', '--force', path, url+'/'+str(path.relative_to(directory)))


def fetch_stage(url, target, names=None):
    target = Path(target)
    marker = load(get(url+'/stage_complete.json', target/'stage_complete.json'))
    if marker['status'] != 'COMPLETE':
        raise ValueError('Upstream stage is incomplete')
    selected = marker['files'] if names is None else {n: marker['files'][n] for n in names}
    for name, record in selected.items():
        if Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError('Unsafe artifact path')
        path = get(url+'/'+name, target/name)
        if path.stat().st_size != record['bytes'] or digest(path) != record['sha256']:
            raise ValueError('Upstream artifact checksum mismatch: '+name)
    return marker


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('pair-train', 'response', 'pair-merge', 'pair-validate', 'pair-close', 'response-merge'))
    parser.add_argument('index', type=int)
    parser.add_argument('system', choices=('OO', 'pp'))
    parser.add_argument('--manifest', type=Path, default=REPO/'campaign_manifest.json')
    parser.add_argument('--scratch', type=Path, required=True)
    args = parser.parse_args()
    manifest = load(args.manifest)
    tasks = [t for t in manifest['tasks'] if t['system'] == args.system]
    base = manifest['outputs'][args.system]
    work = args.scratch.resolve()
    work.mkdir(parents=True, exist_ok=True)
    out = work/'products'
    selected_task = None
    if args.stage in ('pair-train', 'response', 'pair-validate'):
        selected_task = manifest['tasks'][args.index]
        if selected_task['system'] != args.system:
            raise ValueError('DAG system/task mismatch')
        stem = f"{args.system}_{selected_task['file_index']:03d}"
    python = sys.executable
    pairs = REPO/'validation/run3Pairs'
    response = REPO/'validation/run3Unfolding/batch'
    if args.stage in ('pair-train', 'response'):
        input_url = xrd_path(selected_task['input'])
        parsed = urlsplit(input_url)
        stat = subprocess.check_output(['xrdfs', parsed.scheme+'://'+parsed.netloc, 'stat', parsed.path], text=True)
        size = int(re.search(r'^Size:\s+(\d+)', stat, re.M).group(1))
        modified = datetime.strptime(re.search(r'^MTime:\s+(.+)', stat, re.M).group(1).strip(), '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
        identity = selected_task['input_identity']
        if size != identity['size_bytes'] or int(modified.timestamp()) != identity['mtime_ns']//1_000_000_000:
            raise ValueError('Input changed after campaign freeze')
        forest = get(input_url, work/'input.root')
        if args.stage == 'pair-train':
            run(python, pairs/'batch_worker.py', '--manifest', args.manifest, '--task-index', args.index,
                '--output-dir', out, '--scratch', work/'native', '--input-local', forest)
            destination = base+'/pairs/shards/'+stem
        else:
            run(python, response/'response_worker.py', forest, out, '--system', args.system,
                '--namespace', selected_task['sample_namespace'], '--source-input', selected_task['input'],
                '--analysis-repo', REPO, '--scratch-parent', work)
            receipt = load(out/'receipt.json')
            if receipt['skim']['processed_entries'] != selected_task['entry_limit'] or receipt['skim']['input_entries'] != selected_task['entry_limit']:
                raise ValueError('Response does not cover every frozen input entry')
            receipt['frozen_task'] = selected_task
            receipt['campaign_manifest_sha256'] = digest(args.manifest)
            save(out/'receipt.json', receipt)
            destination = base+'/response/shards/'+stem
    elif args.stage == 'pair-merge':
        shards = work/'shards'
        shards.mkdir()
        def collect(task):
            prefix = f"{args.system}_{task['file_index']:03d}"
            local = work/'download'/prefix
            fetch_stage(base+'/pairs/shards/'+prefix, local)
            for path in local.iterdir():
                if path.name.startswith(prefix+'.'):
                    path.rename(shards/path.name)
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(collect, tasks))
        out.mkdir()
        run(python, pairs/'batch_pairs.py', 'merge-training', args.manifest, shards, out/(args.system+'.training.json'),
            '--system', args.system, '--identity-db', work/'identities.sqlite', '--min-pairs', '10')
        destination = base+'/pairs/training'
    elif args.stage == 'pair-validate':
        cache = stem+'.holdout.jsonl.gz'
        fetch_stage(base+'/pairs/shards/'+stem, work/'cache', [cache])
        table = args.system+'.training.json'
        fetch_stage(base+'/pairs/training', work/'training', [table])
        out.mkdir()
        run(python, pairs/'batch_pairs.py', 'validate', work/'cache'/cache, work/'training'/table, out/(stem+'.validation.json.gz'))
        destination = base+'/pairs/validation/'+stem
    elif args.stage == 'pair-close':
        table = args.system+'.training.json'
        fetch_stage(base+'/pairs/training', work/'training', [table])
        validation = work/'validation'
        validation.mkdir()
        def collect(task):
            prefix = f"{args.system}_{task['file_index']:03d}"
            name = prefix+'.validation.json.gz'
            local = work/'download'/prefix
            fetch_stage(base+'/pairs/validation/'+prefix, local, [name])
            (local/name).rename(validation/name)
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(collect, tasks))
        out.mkdir()
        result = out/(args.system+'.pairs.json')
        run(python, pairs/'batch_pairs.py', 'merge-validation', work/'training'/table, validation, result)
        first = f"{args.system}_{tasks[0]['file_index']:03d}"
        cache = first+'.holdout.jsonl.gz'
        fetch_stage(base+'/pairs/shards/'+first, work/'representative', [cache])
        run(python, pairs/'batch_plot.py', result, out/'plots', '--first-then-all', '--events', work/'representative'/cache)
        destination = base+'/pairs/final'
    else:
        dirs = []
        def collect(task):
            prefix = f"{args.system}_{task['file_index']:03d}"
            local = work/'shards'/prefix
            fetch_stage(base+'/response/shards/'+prefix, local)
            receipt = load(local/'receipt.json')
            if receipt.get('frozen_task') != task or receipt.get('campaign_manifest_sha256') != digest(args.manifest):
                raise ValueError('Response shard manifest mismatch')
            return str(local)
        with ThreadPoolExecutor(max_workers=4) as pool:
            dirs = list(pool.map(collect, tasks))
        save(work/'inputs.json', dirs)
        run(python, response/'merge_response.py', work/'inputs.json', out, '--system', args.system, '--manifest', args.manifest)
        destination = base+'/response/final'
    save(out/'campaign_stage.json', {'stage': args.stage, 'system': args.system, 'task_index': args.index,
         'manifest_sha256': digest(args.manifest), 'task': selected_task, 'destination': destination})
    publish(out, destination)
    print(json.dumps({'status': 'CAMPAIGN_STAGE_PUBLISHED', 'stage': args.stage, 'destination': destination}), flush=True)

if __name__ == '__main__':
    main()
