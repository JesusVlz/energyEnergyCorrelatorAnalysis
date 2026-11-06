#!/usr/bin/env python3
"""Snapshot stable OO/ppRef CRAB ROOT outputs using filesystem metadata only.

Select tasks from the current repository manifests, substituting OO
part02/chunk02 v3 for v2. A fresh scan checks each candidate twice and
requires its initial mtime to predate the snapshot start by 120 seconds.
Use --input-inventory to recheck only members of an earlier inventory.
Use --max-oo-part to limit OO parts while retaining every ppRef task.
Ambiguous versions, namespaces, job IDs and retry filenames are not merged.
"""

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import time


DEFAULT_SOURCE = Path('/eos/cms/store/group/phys_heavyions/jvelazqu/eec')
VERSION_OVERRIDE = {'eecAnalysis_OO_mixedPerp_part02_chunk02_Run3_v2':
                    'eecAnalysis_OO_mixedPerp_part02_chunk02_Run3_v3'}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def task_specs(repo):
    manifests = {'OO': repo / 'crab/OOForestfiles/manifest.json',
                 'ppRef': repo / 'crab/ppRefForestfiles/manifest.json'}
    specs = []
    for sample, path in manifests.items():
        manifest = json.loads(path.read_text())
        if manifest['units_per_job'] != 1:
            raise ValueError('Job-ID coverage requires one source file per job')
        for task in manifest['tasks']:
            spec = dict(task, sample=sample, expected_jobs=task['jobs'],
                        manifest_request_name=task['request_name'])
            spec['request_name'] = VERSION_OVERRIDE.get(task['request_name'], task['request_name'])
            if spec['request_name'] != task['request_name']:
                spec['version_override_reason'] = 'OO part02 chunk02 v3 supersedes v2'
            specs.append(spec)
    if len({spec['request_name'] for spec in specs}) != len(specs):
        raise ValueError('Duplicate requests in source manifests')
    return specs, {sample: {'path': str(path), 'sha256': sha256(path)}
                   for sample, path in manifests.items()}


def select_task_scope(specs, max_oo_part=None):
    authoritative_max = max(spec['part'] for spec in specs if spec['sample'] == 'OO')
    if max_oo_part is not None and not 1 <= max_oo_part <= authoritative_max:
        raise ValueError(f'--max-oo-part must be between 1 and {authoritative_max}')
    selected = [spec for spec in specs
                if spec['sample'] != 'OO' or max_oo_part is None
                or spec['part'] <= max_oo_part]
    selected_names = {spec['request_name'] for spec in selected}
    scope = {
        'max_oo_part': max_oo_part,
        'authoritative_max_oo_part': authoritative_max,
        'includes_all_authoritative_tasks': len(selected) == len(specs),
        'selected_task_count': len(selected),
        'authoritative_task_count': len(specs),
        'excluded_authoritative_request_names': sorted(
            spec['request_name'] for spec in specs
            if spec['request_name'] not in selected_names),
    }
    return selected, scope


def validate_candidates(task, records, source):
    root = source / task['request_name']
    pattern = re.compile(re.escape(task['request_name']) + r'_(\d+)\.root$')
    ids = set()
    paths = set()
    namespaces = set()
    for item in records:
        path = Path(item['path'])
        relative = path.relative_to(root)
        match = pattern.fullmatch(path.name)
        if not path.parent.name.isdigit() or not match or len(relative.parts) < 4:
            raise ValueError('Unexpected output or retry filename: ' + str(path))
        job = int(match.group(1))
        if item['job_id'] != job or not 1 <= job <= task['expected_jobs']:
            raise ValueError('Invalid/out-of-range job ID: ' + str(path))
        if item.get('retry') is not None:
            raise ValueError('Explicit retry requires manual selection: ' + str(path))
        if job in ids or str(path) in paths:
            raise ValueError('Duplicate output job ID/path: ' + str(path))
        ids.add(job)
        paths.add(str(path))
        namespaces.add(str(relative.parent.parent))
    if len(namespaces) > 1:
        raise ValueError('Multiple submission namespaces for ' + task['request_name'])


def scan_task(task, source):
    root = source / task['request_name']
    if not root.is_dir():
        raise ValueError('Missing output task directory: ' + str(root))
    records = []
    errors = []
    pattern = re.compile(re.escape(task['request_name']) + r'_(\d+)(?:_(\d+))?\.root$')
    for current, dirs, files in os.walk(root, onerror=lambda error: errors.append(str(error))):
        dirs[:] = sorted(name for name in dirs if name not in ('log', 'logs'))
        for name in sorted(files):
            if not name.endswith('.root'):
                continue
            path = Path(current) / name
            match = pattern.fullmatch(name)
            stat = path.stat()
            relative = path.relative_to(root)
            records.append({'path': str(path),
                            'job_id': int(match.group(1)) if match else None,
                            'retry': int(match.group(2)) if match and match.group(2) else None,
                            'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                            'timestamp_namespace': str(relative.parent.parent),
                            'numeric_directory': path.parent.name})
    if errors:
        raise ValueError('Incomplete directory scan: ' + '; '.join(errors))
    validate_candidates(task, records, source)
    records.sort(key=lambda item: item['job_id'])
    return dict(task, directory=str(root), files=records)


def candidate_tasks(args, specs, manifests):
    if args.input_inventory is None:
        with ThreadPoolExecutor(max_workers=4) as executor:
            tasks = list(executor.map(lambda task: scan_task(task, args.source_root), specs))
        return tasks, None
    report = json.loads(args.input_inventory.read_text())
    if Path(report['base']) != args.source_root:
        raise ValueError('Input inventory source root differs from --source-root')
    for sample, info in manifests.items():
        if report['source_manifests'][sample]['sha256'] != info['sha256']:
            raise ValueError('Source manifest changed since the input inventory: ' + sample)
    found = {task['request_name']: task for task in report['tasks']}
    if len(found) != len(report['tasks']) or set(found) != {spec['request_name'] for spec in specs}:
        raise ValueError('Input inventory tasks differ from the selected authoritative tasks; '
                         'use the same --max-oo-part selection as the input inventory')
    tasks = []
    for spec in specs:
        old = found[spec['request_name']]
        if old['expected_jobs'] != spec['expected_jobs'] or old.get('metadata_errors'):
            raise ValueError('Incomplete/incompatible input inventory: ' + spec['request_name'])
        for key in ('duplicate_job_ids', 'explicit_retry_files', 'unrecognized_names',
                    'root_files_outside_numeric_directories', 'out_of_range_job_ids'):
            if old.get(key):
                raise ValueError('Input inventory contains ambiguous outputs: ' + key)
        validate_candidates(spec, old['files'], args.source_root)
        tasks.append(dict(spec, directory=str(args.source_root / spec['request_name']),
                          files=old['files']))
    return tasks, {'path': str(args.input_inventory), 'sha256': sha256(args.input_inventory),
                   'original_summary': report.get('summary'),
                   'original_scan_started_utc': report.get('scan_started_utc'),
                   'original_scan_finished_utc': report.get('scan_finished_utc')}


def stabilize(task, cutoff_ns):
    original = task['files']
    ready = []
    withheld = []
    for item in original:
        reasons = []
        check = None
        try:
            stat = Path(item['path']).stat()
            check = {'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
            if item['bytes'] <= 0 or stat.st_size <= 0:
                reasons.append('zero_size_at_one_or_both_checks')
            if stat.st_size != item['bytes'] or stat.st_mtime_ns != item['mtime_ns']:
                reasons.append('size_or_mtime_changed_between_checks')
            if item['mtime_ns'] > cutoff_ns:
                reasons.append('younger_than_minimum_age_at_snapshot_start')
        except OSError as error:
            reasons.append('restat_failed: ' + str(error))
        if reasons:
            withheld.append(dict(item, reasons=reasons, recheck=check))
        else:
            ready.append(item)
    ids = {item['job_id'] for item in ready}
    original_ids = {item['job_id'] for item in original}
    namespaces = defaultdict(lambda: {'files': 0, 'bytes': 0, 'numeric_directories': set()})
    for item in ready:
        namespace = namespaces[item['timestamp_namespace']]
        namespace['files'] += 1
        namespace['bytes'] += item['bytes']
        namespace['numeric_directories'].add(item['numeric_directory'])
    for namespace in namespaces.values():
        namespace['numeric_directories'] = sorted(namespace['numeric_directories'])
    result = {key: value for key, value in task.items() if key != 'files'}
    result.update({'directory_exists': True,
                   'original_observed_files': len(original),
                   'original_total_bytes': sum(item['bytes'] for item in original),
                   'original_unique_job_ids': len(original_ids),
                   'original_missing_jobs': task['expected_jobs'] - len(original_ids),
                   'observed_files': len(ready), 'unique_job_ids': len(ids),
                   'total_bytes': sum(item['bytes'] for item in ready),
                   'missing_job_ids': sorted(set(range(1, task['expected_jobs'] + 1)) - ids),
                   'timestamp_namespaces': dict(sorted(namespaces.items())),
                   'withheld_files': withheld, 'withheld_file_count': len(withheld),
                   'duplicate_job_ids': {}, 'explicit_retry_files': [],
                   'unrecognized_names': [], 'zero_byte_files': [],
                   'root_files_outside_numeric_directories': [],
                   'out_of_range_job_ids': [], 'metadata_errors': [], 'files': ready})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New inventory JSON path')
    parser.add_argument('--input-inventory', type=Path,
                        help='Recheck only existing members of this earlier inventory')
    parser.add_argument('--source-root', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--max-oo-part', type=int,
                        help='Include OO parts up to this number and all ppRef tasks (default: all)')
    parser.add_argument('--min-age-seconds', type=float, default=120,
                        help='Minimum mtime age at snapshot start (default: 120)')
    args = parser.parse_args()
    if args.min_age_seconds < 0:
        parser.error('--min-age-seconds must be nonnegative')
    if args.output.exists():
        parser.error('--output must not already exist; use a new snapshot filename')
    start_ns = time.time_ns()
    started = utc_now()
    authoritative_specs, manifests = task_specs(args.repo_root)
    try:
        specs, selection_scope = select_task_scope(authoritative_specs, args.max_oo_part)
    except ValueError as error:
        parser.error(str(error))
    full_campaign_expected_jobs = {
        sample: sum(spec['expected_jobs'] for spec in authoritative_specs if spec['sample'] == sample)
        for sample in ('OO', 'ppRef')
    }
    candidates, previous = candidate_tasks(args, specs, manifests)
    cutoff_ns = start_ns - int(args.min_age_seconds * 1_000_000_000)
    with ThreadPoolExecutor(max_workers=4) as executor:
        tasks = list(executor.map(lambda task: stabilize(task, cutoff_ns), candidates))
    summary = {}
    for sample in ('OO', 'ppRef'):
        group = [task for task in tasks if task['sample'] == sample]
        summary[sample] = {
            'tasks': len(group), 'expected_jobs': sum(task['expected_jobs'] for task in group),
            'full_campaign_expected_jobs': full_campaign_expected_jobs[sample],
            'excluded_authoritative_jobs': (full_campaign_expected_jobs[sample]
                                            - sum(task['expected_jobs'] for task in group)),
            'original_observed_files': sum(task['original_observed_files'] for task in group),
            'original_total_bytes': sum(task['original_total_bytes'] for task in group),
            'observed_files': sum(task['observed_files'] for task in group),
            'unique_task_job_ids': sum(task['unique_job_ids'] for task in group),
            'missing_jobs': sum(len(task['missing_job_ids']) for task in group),
            'total_bytes': sum(task['total_bytes'] for task in group),
            'withheld_files': sum(task['withheld_file_count'] for task in group),
            'duplicate_task_job_ids': 0, 'metadata_errors': 0,
        }
    selected = {task['request_name'] for task in tasks}
    known = sorted(path.name for path in args.source_root.glob('eecAnalysis_*mixedPerp*') if path.is_dir())
    report = {
        'snapshot_started_utc': started, 'snapshot_finished_utc': utc_now(),
        'snapshot_started_ns': start_ns, 'min_age_seconds': args.min_age_seconds,
        'base': str(args.source_root), 'metadata_only': True, 'stable_snapshot': True,
        'selection_policy': 'Manifest tasks with explicit OO part02/chunk02 v3 override; optional OO part limit retains all ppRef tasks; exclude unsplit v1 and other unselected tasks',
        'selection_scope': selection_scope,
        'full_campaign_expected_jobs': full_campaign_expected_jobs,
        'stability_policy': 'Keep only nonzero files with unchanged size and mtime at both checks, older than minimum age at snapshot start; input inventory mode never adds members',
        'source_manifests': manifests, 'input_inventory': previous,
        'observed_task_directories': known,
        'excluded_task_directories': [name for name in known if name not in selected],
        'summary': summary, 'tasks': tasks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    for sample, values in summary.items():
        print(f"{sample}: {values['observed_files']:,}/{values['expected_jobs']:,} ready files; "
              f"{values['withheld_files']:,} withheld; {values['total_bytes']:,} bytes")
    print('Inventory: ' + str(args.output))


if __name__ == '__main__':
    main()
