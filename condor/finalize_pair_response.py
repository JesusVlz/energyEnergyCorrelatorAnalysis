#!/usr/bin/env python3
"""Verify every combined forest archive, then run the existing finalizers.

This helper does not submit jobs or publish outputs. Pair training counts are
summed globally before normalization; validation uses that single global table.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile

REPO = Path(__file__).resolve().parents[1]
VERSION = 'pair-response-finalizer-v1'


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def relative_path(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or '..' in path.parts or '\\' in name or not path.parts:
        raise ValueError('Unsafe archive/artifact path: ' + repr(name))
    return path


def prefix(task):
    return '{}_{:03d}'.format(task['system'], task['file_index'])


def expected_tasks(manifest, system):
    tasks = [task for task in manifest['tasks'] if task['system'] == system]
    if not tasks:
        raise ValueError('No expected tasks for ' + system)
    if len({prefix(task) for task in tasks}) != len(tasks):
        raise ValueError('Duplicate task prefix in frozen manifest')
    inputs = [str(task['input']).replace('/eos/home-j/', '/eos/user/j/') for task in tasks]
    if len(set(inputs)) != len(tasks):
        raise ValueError('Duplicate source input in frozen manifest')
    for task in tasks:
        if task['entry_limit'] <= 0 or task['input_identity']['entries'] != task['entry_limit']:
            raise ValueError('Task does not cover a complete nonempty forest')
    if 'counts' in manifest and manifest['counts'][system] != len(tasks):
        raise ValueError('Frozen system task count mismatch')
    return tasks


def extract_checked(archive_path, target, kind, task, manifest_sha256):
    """Extract only one kind; never allow tar links, traversal or hidden files."""
    target = Path(target)
    if target.exists():
        raise ValueError('Archive extraction target must be new')
    with tarfile.open(archive_path, 'r:gz') as archive:
        names, selected = set(), {}
        for member in archive.getmembers():
            path = relative_path(member.name)
            if path.parts[0] not in ('pairs', 'response'):
                raise ValueError('Unexpected archive root: ' + member.name)
            normalized = str(path)
            if normalized in names:
                raise ValueError('Duplicate archive member: ' + normalized)
            names.add(normalized)
            if not member.isdir() and not member.isfile():
                raise ValueError('Archive links/special files are forbidden: ' + normalized)
            if path.parts[0] == kind and member.isfile():
                if len(path.parts) < 2:
                    raise ValueError('Archive kind must be a directory')
                selected[str(PurePosixPath(*path.parts[1:]))] = member
        if 'stage_complete.json' not in selected:
            raise ValueError('Missing selected-kind completion marker')
        marker = json.load(archive.extractfile(selected['stage_complete.json']))
        if (marker.get('status') != 'COMPLETE' or marker.get('frozen_task') != task or
                marker.get('campaign_manifest_sha256') != manifest_sha256):
            raise ValueError('Incomplete or mismatched frozen archive task/manifest')
        files = marker['files']
        if set(selected) != set(files) | {'stage_complete.json'}:
            raise ValueError('Archive files differ from completion marker')
        for name, record in files.items():
            relative_path(name)
            if (not isinstance(record['bytes'], int) or record['bytes'] < 0 or
                    not re.fullmatch(r'[a-f0-9]{64}', record['sha256']) or
                    selected[name].size != record['bytes']):
                raise ValueError('Invalid artifact size/checksum record: ' + name)
        target.mkdir(parents=True)
        for name, member in selected.items():
            path = target / name
            path.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(member) as source, path.open('xb') as destination:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
            if name in files and (path.stat().st_size != files[name]['bytes'] or
                                  digest(path) != files[name]['sha256']):
                raise ValueError('Artifact checksum mismatch: ' + name)
    provenance = load(target / 'combined_worker_provenance.json')
    audit = provenance.get('shared_identity_audit', {})
    if (provenance.get('status') != 'COMPLETE' or provenance.get('frozen_task') != task or
            provenance.get('campaign_manifest_sha256') != manifest_sha256 or
            audit.get('status') != 'PASS' or audit.get('matched_raw_physical_ids') != task['entry_limit']):
        raise ValueError('Combined worker provenance/identity audit mismatch')
    receipt = load(target / (prefix(task) + '.receipt.json' if kind == 'pairs' else 'receipt.json'))
    if kind == 'pairs':
        if receipt.get('task') != task or receipt['event_counts']['raw'] != task['entry_limit']:
            raise ValueError('Pair receipt task/coverage mismatch')
    else:
        if (receipt.get('status') != 'COMPLETE_ENGINEERING_ONLY' or receipt.get('frozen_task') != task or
                receipt.get('campaign_manifest_sha256') != manifest_sha256 or
                receipt['skim']['input_entries'] != task['entry_limit'] or
                receipt['skim']['processed_entries'] != task['entry_limit']):
            raise ValueError('Response receipt task/coverage mismatch')
    return provenance


def fetch(url, destination):
    subprocess.run(['xrdcp', '--retry', '3', '--cksum', 'adler32', url, str(destination)], check=True)


def collect_inputs(manifest, manifest_sha256, system, kind, work):
    tasks = expected_tasks(manifest, system)
    base = manifest['outputs'][system].rstrip('/')
    if not base.startswith('root://'):
        raise ValueError('Expected an explicit XRootD output base')
    records, folders, signature = [], [], None
    for task in tasks:
        name = prefix(task)
        archive = work / (name + '.tar.gz')
        url = base + '/PairResponse/' + name + '.tar.gz'
        fetch(url, archive)
        folder = work / 'shards' / name
        provenance = extract_checked(archive, folder, kind, task, manifest_sha256)
        actual = {'source_fingerprints': provenance['source_fingerprints'],
                  'driver_sha256': provenance['driver_sha256']}
        if signature is not None and actual != signature:
            raise ValueError('Mixed combined-worker source/binary/card versions')
        signature = actual
        records.append({'prefix': name, 'url': url, 'archive_bytes': archive.stat().st_size,
                        'archive_sha256': digest(archive),
                        'completion_sha256': digest(folder / 'stage_complete.json'),
                        'raw_entries': task['entry_limit']})
        folders.append(folder)
        archive.unlink()  # Keep only the requested kind in local scratch.
    report = {'status': 'ALL_EXPECTED_ARCHIVES_VERIFIED', 'system': system, 'kind': kind,
              'manifest_sha256': manifest_sha256, 'files': len(tasks),
              'raw_entries': sum(task['entry_limit'] for task in tasks),
              'source_signature': signature, 'archives': records}
    save(work / 'verified_inputs.json', report)
    return tasks, folders, report


def run_pipeline(manifest_path, manifest, system, kind, tasks, folders, work, output):
    commands = []

    def run(*args):
        command = [str(arg) for arg in args]
        commands.append(command)
        subprocess.run(command, cwd=REPO, check=True)

    python = sys.executable
    if kind == 'pairs':
        output.mkdir(parents=True)
        shards, validation = work / 'pair_shards', work / 'validation'
        shards.mkdir()
        validation.mkdir()
        for task, folder in zip(tasks, folders):
            for path in folder.iterdir():
                if path.is_file() and path.name.startswith(prefix(task) + '.'):
                    os.link(path, shards / path.name)
        helper = REPO / 'validation/run3Pairs/batch_pairs.py'
        training, result = output / (system + '.training.json'), output / (system + '.pairs.json')
        run(python, helper, 'merge-training', manifest_path, shards, training, '--system', system,
            '--identity-db', work / 'identities.sqlite', '--min-pairs', manifest.get('scope', {}).get('pair_min_pairs', 10))
        for task in tasks:
            name = prefix(task)
            run(python, helper, 'validate', shards / (name + '.holdout.jsonl.gz'), training,
                validation / (name + '.validation.json.gz'))
        run(python, helper, 'merge-validation', training, validation, result)
        run(python, REPO / 'validation/run3Pairs/batch_plot.py', result, output / 'plots',
            '--first-then-all', '--events', shards / (prefix(tasks[0]) + '.holdout.jsonl.gz'))
        if load(result).get('batch_status') != 'ALL_EXPECTED_SHARDS_TRAINED_AND_HELD_OUT':
            raise ValueError('Pair finalizer did not complete all expected shards')
    else:
        inputs = work / 'response_inputs.json'
        save(inputs, [str(folder) for folder in folders])
        run(python, REPO / 'validation/run3Unfolding/batch/merge_response.py', inputs, output,
            '--system', system, '--manifest', manifest_path, '--scratch-parent', work)
        if load(output / 'merge_receipt.json').get('status') != 'COMPLETE_RESPONSE_MATRICES_ENGINEERING_ONLY':
            raise ValueError('Response finalizer did not complete')
    return commands


def finalize(args):
    manifest_path = args.manifest.resolve(strict=True)
    manifest = load(manifest_path)
    work, output = args.work_dir.resolve(), args.output_dir.resolve()
    if output.exists():
        raise ValueError('Final output directory must be new')
    if work.exists() and any(work.iterdir()):
        raise ValueError('Finalizer work directory must be new or empty')
    if str(work).startswith(('/eos/', '/afs/')):
        raise ValueError('Finalizer work directory must use worker-local scratch')
    work.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    tasks, folders, report = collect_inputs(manifest, digest(manifest_path), args.system, args.kind, work)
    commands = run_pipeline(manifest_path, manifest, args.system, args.kind, tasks, folders, work, output)
    helper_paths = (['validation/run3Pairs/batch_pairs.py', 'validation/run3Pairs/pair_counts.py',
                     'validation/run3Pairs/batch_plot.py', 'validation/run3Pairs/plot_closure.py']
                    if args.kind == 'pairs' else
                    ['validation/run3Unfolding/batch/merge_response.py',
                     'validation/run3Unfolding/batch/response_worker.py',
                     'validation/run3Unfolding/project_response.py'])
    report.update(status='COMPLETE_PAIR_RESPONSE_FINALIZATION', commands=commands,
                  finalizer_version=VERSION, finalizer_sha256=digest(__file__),
                  postprocessing_sha256={name: digest(REPO / name) for name in helper_paths},
                  scope='Existing raw-count merges, holdout closure or response projection; no unfolding',
                  outputs={str(path.relative_to(output)): {'bytes': path.stat().st_size, 'sha256': digest(path)}
                           for path in sorted(output.rglob('*')) if path.is_file()})
    save(output / 'finalization_receipt.json', report)
    print(json.dumps({'status': report['status'], 'system': args.system, 'kind': args.kind,
                      'files': report['files'], 'raw_entries': report['raw_entries']}), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--system', choices=('OO', 'pp'), required=True)
    parser.add_argument('--kind', choices=('pairs', 'response'), required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    finalize(parser.parse_args())


if __name__ == '__main__':
    main()
