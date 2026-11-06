#!/usr/bin/env python3
"""Prepare a frozen OO/ppRef CRAB-data merge DAG; never submit jobs.

Inventory comes from inventoryDataOutputs.py. Raw files are merged in bounded
batches, batches into their CRAB parts, then parts into separate sample totals.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import urlsplit

CMSSW = '/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_0_11'
SAMPLES = {'OO': 9, 'ppRef': 11}
HELPERS = ('merge_data_worker.py', 'merge_root.py', 'frozen_input.py', 'prepareDataMerge.py', 'inventoryDataOutputs.py')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    checksum = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def safe_path(path):
    value = str(Path(path).absolute())
    require(re.fullmatch(r'/[A-Za-z0-9_./-]+', value), 'Unsafe submission path: ' + value)
    return Path(value)


def eos_url(path):
    path = str(path)
    require(path.startswith('/eos/cms/store/'), 'Expected CMS EOS storage path')
    return 'root://eoscms.cern.ch/' + path


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def disk_mb(input_bytes):
    return max(16000, math.ceil(3.9 * input_bytes / 1000000) + 2000)


def resources(node, input_bytes=None):
    return {'cpus': 1, 'memory_mb': {'batch': 8000, 'part': 16000, 'sample': 24000}[node['kind']],
            'disk_mb': disk_mb(input_bytes) if input_bytes is not None else 64000,
            'job_flavour': 'workday' if node['kind'] == 'batch' else 'tomorrow'}


def write_resources(path, record):
    path.write_text(f"request_cpus = {record['cpus']}\nrequest_memory = {record['memory_mb']}MB\n"
                    f"request_disk = {record['disk_mb']}MB\n+JobFlavour = \"{record['job_flavour']}\"\n")


def update_resources(node_path):
    node = json.loads(node_path.read_text())
    require(node['kind'] in ('part', 'sample'), 'Only reductions need dynamic scratch sizing')
    require(node['sample'] in SAMPLES and node['data_type'] == SAMPLES[node['sample']],
            'Invalid reduction sample/DataType')
    parents = node.get('parents')
    require(isinstance(parents, list) and bool(parents) and
            len(parents) == len(node['inputs']), 'Parent/input count mismatch')
    require(all(isinstance(parent_id, str) and
                re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', parent_id) for parent_id in parents) and
            len(parents) == len(set(parents)), 'Invalid or duplicate parent node IDs')
    workspace = node_path.parent.parent
    total = 0
    for parent_id, entry in zip(parents, node['inputs']):
        parsed = urlsplit(entry['report_url'])
        require(parsed.scheme == 'root' and parsed.netloc == 'eoscms.cern.ch', 'Unexpected report host')
        require(('/' + parsed.path.lstrip('/')).startswith('/eos/cms/store/'), 'Unexpected report path')
        # DAG PRE scripts run on the scheduler. Condor already returned these
        # reports to AFS before marking the parent jobs successful, so sizing
        # must not depend on an EOS mount or credentials on that host.
        parent = json.loads((workspace / 'reports' / (parent_id + '.json')).read_text())
        require(parent.get('node_id') == parent_id, 'Wrong parent node ID: ' + parent_id)
        require(parent.get('status') in ('PASS', 'PARTIAL') and
                parent.get('sample') == node['sample'] and
                parent.get('data_type') == node['data_type'], 'Invalid parent report')
        publication = parent.get('publication', {})
        require(publication.get('root_url') == entry['url'] and
                publication.get('report_url') == entry['report_url'] and
                publication.get('root_published') is True and
                publication.get('report_published') is True, 'Wrong/unpublished parent')
        size = parent['artifact']['size_bytes']
        require(type(size) is int and size > 0, 'Invalid parent artifact size')
        total += size
    record = resources(node, total)
    target = workspace / 'resources' / (node['node_id'] + '.sub')
    temporary = target.with_suffix('.tmp')
    write_resources(temporary, record)
    temporary.replace(target)
    print(json.dumps({'node': node['node_id'], 'parent_bytes': total, 'resources': record}))


def build_plan(inventory, output_dir, batch_size=500):
    require(batch_size > 0, 'Batch size must be positive')
    require(inventory.get('stable_snapshot') is True, 'Inventory must contain a stable ready_snapshot')
    nodes, part_nodes, seen, task_ids = [], {sample: [] for sample in SAMPLES}, set(), set()
    for task in inventory['tasks']:
        sample = task['sample']
        require(sample in SAMPLES, 'Unexpected sample')
        task_id = f"{sample}_part{task['part']:02d}_chunk{task['chunk']:02d}"
        require(task_id not in task_ids, 'Repeated task partition')
        task_ids.add(task_id)
        for key in ('duplicate_job_ids', 'explicit_retry_files', 'unrecognized_names',
                    'out_of_range_job_ids', 'metadata_errors', 'root_files_outside_numeric_directories'):
            require(not task.get(key), 'Inventory anomaly in ' + task_id + ': ' + key)
        rows = task['files']
        require(rows, 'No stable files for ' + task_id)
        ids = [row['job_id'] for row in rows]
        require(len(ids) == len(set(ids)), 'Duplicate job IDs in ' + task_id)
        require(set(ids) <= set(range(1, task['expected_jobs'] + 1)), 'Out-of-range job IDs')
        require(task['missing_job_ids'] == sorted(set(range(1, task['expected_jobs'] + 1)) - set(ids)), 'Missing job IDs disagree with snapshot files')
        for row in rows:
            require(row['bytes'] > 0 and row['path'].startswith(task['directory'].rstrip('/') + '/'), 'Invalid raw input')
            require(row['path'] not in seen, 'Repeated physical input path')
            seen.add(row['path'])
        snapshot = {'request_name': task['request_name'], 'expected_jobs': task['expected_jobs'],
                    'observed_jobs': len(rows), 'missing_job_ids': task['missing_job_ids'],
                    'complete': len(rows) == task['expected_jobs'],
                    'scan_started_utc': inventory['snapshot_started_utc'],
                    'withheld_files': task.get('withheld_files', [])}
        batches = []
        for number, start in enumerate(range(0, len(rows), batch_size), 1):
            entries = rows[start:start + batch_size]
            node_id = f'{task_id}_batch{number:03d}'
            node = {'schema_version': 1, 'node_id': node_id, 'sample': sample, 'data_type': SAMPLES[sample],
                    'kind': 'batch', 'inputs': [dict(url=eos_url(row['path']), size_bytes=row['bytes'],
                        job_id=row['job_id'], request_name=task['request_name'], mtime_ns=row['mtime_ns']) for row in entries],
                    'output_url': eos_url(output_dir / sample / 'batches' / (node_id + '.root')),
                    'report_url': eos_url(output_dir / sample / 'batches' / (node_id + '.json')),
                    'parents': [], 'snapshot': {**snapshot, 'batch_job_ids': [r['job_id'] for r in entries]}}
            node['resources'] = resources(node, sum(row['bytes'] for row in entries))
            nodes.append(node)
            batches.append(node)
        part = {'schema_version': 1, 'node_id': task_id, 'sample': sample, 'data_type': SAMPLES[sample],
                'kind': 'part', 'inputs': [dict(url=n['output_url'], report_url=n['report_url']) for n in batches],
                'output_url': eos_url(output_dir / sample / 'parts' / (task_id + '.root')),
                'report_url': eos_url(output_dir / sample / 'parts' / (task_id + '.json')),
                'parents': [n['node_id'] for n in batches], 'snapshot': snapshot}
        part['resources'] = resources(part)
        nodes.append(part)
        part_nodes[sample].append(part)
    for sample, parents in part_nodes.items():
        require(parents, 'No parts for ' + sample)
        node = {'schema_version': 1, 'node_id': 'merged_' + sample, 'sample': sample, 'data_type': SAMPLES[sample],
                'kind': 'sample', 'inputs': [dict(url=p['output_url'], report_url=p['report_url']) for p in parents],
                'output_url': eos_url(output_dir / ('merged_' + sample + '.root')),
                'report_url': eos_url(output_dir / ('merged_' + sample + '.json')),
                'parents': [p['node_id'] for p in parents],
                'snapshot': {'scan_started_utc': inventory['snapshot_started_utc'],
                             'complete': all(p['snapshot']['complete'] for p in parents),
                             'expected_jobs': sum(p['snapshot']['expected_jobs'] for p in parents),
                             'observed_jobs': sum(p['snapshot']['observed_jobs'] for p in parents),
                             'parts': [p['snapshot'] for p in parents]}}
        node['snapshot']['full_campaign_expected_jobs'] = inventory.get(
            'full_campaign_expected_jobs', {}).get(sample, inventory['summary'][sample]['expected_jobs'])
        if 'selection_scope' in inventory:
            node['snapshot']['selection_scope'] = inventory['selection_scope']
        node['resources'] = resources(node)
        nodes.append(node)
    require(len({n['node_id'] for n in nodes}) == len(nodes), 'Duplicate node ID')
    return nodes


def prepare(inventory_path, workspace, output_dir, tag, batch_size=500):
    workspace, output_dir = safe_path(workspace), safe_path(output_dir)
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', tag), 'Invalid tag')
    require(not workspace.exists(), 'Workspace already exists: ' + str(workspace))
    require(not output_dir.exists(), 'Output snapshot already exists: ' + str(output_dir))
    inventory = json.loads(inventory_path.read_text())
    nodes = build_plan(inventory, output_dir, batch_size)
    source = Path(__file__).resolve().parent
    for helper in HELPERS:
        require((source / helper).is_file(), 'Missing helper ' + helper)
    workspace.mkdir(parents=True)
    for name in ('nodes', 'jobs', 'resources', 'logs', 'reports', 'checks'):
        (workspace / name).mkdir()
    for helper in HELPERS:
        shutil.copy2(source / helper, workspace / helper)
    shutil.copy2(inventory_path, workspace / 'inventory.json')
    helper_hashes = {name: sha256(workspace / name) for name in HELPERS}
    wrapper = '''#!/bin/bash
set -Eeuo pipefail
scratch=$(pwd -P)
case "$scratch" in /afs/*|/eos/*) echo 'Local scratch required' >&2; exit 2 ;; esac
printf '{"status":"FAIL","error":"worker did not start"}\\n' > result.json
: > good_files.txt
: > bad_files.txt
[[ $# == 2 && $1 =~ ^[A-Za-z0-9_.-]+[.]json$ && $2 =~ ^[a-f0-9]{64}$ ]]
printf '%s  %s\\n' "$2" "$1" | sha256sum -c -
'''
    for name in ('merge_data_worker.py', 'merge_root.py', 'frozen_input.py'):
        wrapper += f"echo '{helper_hashes[name]}  {name}' | sha256sum -c -\n"
    wrapper += f'''set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cd {CMSSW}
runtime=$(scram runtime -sh)
eval "$runtime"
set -u
cd "$scratch"
set +e
python3 merge_data_worker.py --manifest "$1" --work-dir "$scratch/merge-work"
status=$?
set -e
for name in result.json good_files.txt bad_files.txt; do
  if [[ -f merge-work/$name ]]; then cp "merge-work/$name" "$name"; fi
done
exit "$status"
'''
    (workspace / 'run_merge.sh').write_text(wrapper)
    (workspace / 'run_merge.sh').chmod(0o755)
    dag = []
    node_records = []
    for node in nodes:
        node_id = node['node_id']
        node_path = workspace / 'nodes' / (node_id + '.json')
        write_json(node_path, node)
        node_hash = sha256(node_path)
        write_resources(workspace / 'resources' / (node_id + '.sub'), node['resources'])
        transfers = ','.join(str(workspace / name) for name in ('merge_data_worker.py', 'merge_root.py', 'frozen_input.py')) + ',' + str(node_path)
        submit = f'''universe = vanilla
initialdir = {workspace}
executable = {workspace}/run_merge.sh
arguments = {node_id}.json {node_hash}
log = {workspace}/logs/jobs.$(ClusterId).log
output = {workspace}/logs/{node_id}.$(ClusterId).$(ProcId).out
error = {workspace}/logs/{node_id}.$(ClusterId).$(ProcId).err
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {transfers}
transfer_output_files = result.json,good_files.txt,bad_files.txt
transfer_output_remaps = "result.json=reports/{node_id}.json;good_files.txt=reports/{node_id}.good.txt;bad_files.txt=reports/{node_id}.bad.txt"
include : {workspace}/resources/{node_id}.sub
requirements = (OpSysAndVer =?= "AlmaLinux9") && (Arch =?= "X86_64")
+JobBatchName = "Run3Data_merge_{tag}"
queue 1
'''
        (workspace / 'jobs' / (node_id + '.sub')).write_text(submit)
        dag.append(f'JOB {node_id} jobs/{node_id}.sub')
        if node['parents']:
            dag.append(f'SCRIPT PRE {node_id} /usr/bin/python3 {workspace}/prepareDataMerge.py --update-resources {node_path}')
            dag.append('PARENT ' + ' '.join(node['parents']) + ' CHILD ' + node_id)
        node_records.append({key: node[key] for key in ('node_id', 'sample', 'kind', 'parents', 'output_url', 'report_url', 'resources')})
        node_records[-1].update(manifest_sha256=node_hash, immediate_inputs=len(node['inputs']))
    (workspace / 'merge.dag').write_text('\n'.join(dag) + '\n')
    submit_sh = '''#!/bin/bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ ${1:-} == --dry-run ]]; then
  [[ $# == 1 ]] || exit 2
  condor_submit_dag -no_submit -MaxJobs 20 merge.dag
  for file in jobs/*.sub; do
    condor_submit -dry-run "checks/$(basename "$file" .sub).ads" "$file"
  done
elif [[ $# == 0 ]]; then
  exec condor_submit_dag -update_submit -MaxJobs 20 merge.dag
else
  echo "Usage: $0 [--dry-run]" >&2
  exit 2
fi
'''
    (workspace / 'submit.sh').write_text(submit_sh)
    (workspace / 'submit.sh').chmod(0o755)
    summary = {'tag': tag, 'workspace': str(workspace), 'output_dir': str(output_dir),
               'snapshot_summary': inventory['summary'], 'inventory_sha256': sha256(inventory_path),
               'helper_sha256': helper_hashes, 'batch_size': batch_size, 'node_count': len(nodes),
               'node_kinds': dict(Counter(n['kind'] for n in nodes)), 'nodes': node_records,
               'submitted': False,
               'full_campaign_expected_jobs': {n['sample']: n['snapshot']['full_campaign_expected_jobs']
                                               for n in nodes if n['kind'] == 'sample'},
               'snapshot_is_complete_campaign': (inventory.get('selection_scope', {}).get(
                   'includes_all_authoritative_tasks', True) and
                   all(n['snapshot']['complete'] for n in nodes if n['kind'] == 'sample'))}
    if 'selection_scope' in inventory:
        summary['selection_scope'] = inventory['selection_scope']
    write_json(workspace / 'manifest.json', summary)
    lines = [f'OO and ppRef data merge snapshot: {tag}', '', 'Prepared only; no jobs submitted.',
             f'Source scan: {inventory["snapshot_started_utc"]}', f'Output: {output_dir}', '',
             'Frozen ready inputs (later CRAB outputs are not added automatically):']
    for sample in SAMPLES:
        end = next(n for n in nodes if n['node_id'] == 'merged_' + sample)
        lines.append(f'  {sample}: {end["snapshot"]["observed_jobs"]:,} / {end["snapshot"]["expected_jobs"]:,} selected expected jobs; '
                     f'{end["snapshot"]["full_campaign_expected_jobs"]:,} full-campaign expected jobs')
    scope = inventory.get('selection_scope', {})
    if scope.get('max_oo_part') is not None:
        lines.append(f'Requested OO coverage: parts 1-{scope["max_oo_part"]} '
                     f'of {scope["authoritative_max_oo_part"]} authoritative parts; ppRef coverage unchanged.')
    lines += ['', f'DAG: {summary["node_kinds"]}; at most 20 Condor nodes submitted simultaneously.',
              'Each raw batch contains at most 500 files (or --batch-size at preparation).',
              'Batch -> CRAB part -> separate OO/ppRef totals; missing job IDs and withheld files are in inventory.json.',
              'OO part02/chunk02 uses v3; superseded v1 and duplicate campaign generations are excluded.', '',
              'Run from this directory:', '  ./submit.sh --dry-run   # parse all submit files; no queue submission',
              '  ./submit.sh             # submit the dependency DAG',
              '  condor_q -dag', '',
              'Workers use CVMFS CMSSW_15_0_11 and local scratch; ROOT files are staged/published using XRootD.',
              'Requests: one CPU; batch 8 GB/workday, part 16 GB/tomorrow, sample 24 GB/tomorrow.',
              'Scratch starts at >=16 GB for batches. A DAG PRE script sizes each reduction from its actual parent artifact sizes.',
              'CERN batch credentials are used. No copied login proxy, AFS CMSSW workarea, or EOS executable is required.', '',
              'Like merge4.sh, raw zero-byte/unreadable/zombie/recovered ROOT files are excluded and explicitly recorded.',
              'Transport/size errors, incompatible cards/axes, and histogram validation failures stop the job.',
              'ROOT histogram bins/errors and exposures are checked; JCard is restored once after hadd.',
              'Successful parent reports and SHA256 are checked, and original source overlap is rejected.',
              'Per-node good/bad lists and result reports return to reports/ even on worker failure.',
              'Successful ROOT and JSON files are published to EOS, with JSON last. Existing outputs are never overwritten.',
              'No automatic retries: inspect failed nodes and publication state before resubmitting a rescue DAG.', '',
              'Final products (consume only with the corresponding successful JSON report):',
              f'  {output_dir}/merged_OO.root', f'  {output_dir}/merged_OO.json',
              f'  {output_dir}/merged_ppRef.root', f'  {output_dir}/merged_ppRef.json', '',
              'These are partial campaign snapshots, not full-production totals. Root integrity does not establish physics closure.',
              'For a later snapshot, rerun inventoryDataOutputs.py with a new output filename and prepareDataMerge.py with a new tag/directory.']
    (workspace / 'README.txt').write_text('\n'.join(lines) + '\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('nodes', 'helper_sha256')}, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path)
    parser.add_argument('--submit-dir', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--tag')
    parser.add_argument('--batch-size', type=int, default=500)
    parser.add_argument('--update-resources', type=Path)
    args = parser.parse_args()
    if args.update_resources:
        update_resources(args.update_resources)
    else:
        require(all((args.inventory, args.submit_dir, args.output_dir, args.tag)), 'Specify --inventory, --submit-dir, --output-dir and --tag')
        prepare(args.inventory, args.submit_dir, args.output_dir, args.tag, args.batch_size)


if __name__ == '__main__':
    main()
