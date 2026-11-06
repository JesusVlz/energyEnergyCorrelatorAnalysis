#!/usr/bin/env python3
"""Reserve a ppRef mixing pool and pack maximum-size CRAB tasks.

Reserve 35 evenly spaced source-list midpoints for mixing, then pack the
remaining files in source order into tasks of 10,000 one-file jobs, with a
smaller final task if needed. Existing identical outputs are checked;
changed outputs are never overwritten.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit


POOL_FILES = 35
MAX_JOBS_PER_TASK = 10000
REDIRECTOR = 'root://xrootd-vanderbilt.sites.opensciencegrid.org/'
LIMIT_SOURCE = (
    'https://github.com/dmwm/CRABServer/blob/'
    '8069e0cb72445e9a3a13dfb98469a0920210ed12/'
    'src/python/TaskWorker/Actions/Splitter.py#L23'
)


def read_list(path):
    return [line.strip() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def logical_file_name(value):
    path = urlsplit(value).path if '://' in value else value
    offset = path.find('/store/')
    if offset < 0:
        raise ValueError('Expected a /store/ file name: ' + value)
    return path[offset:]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def dataset_name(lfn):
    components = lfn.split('/')
    try:
        return components[components.index('HiForestpp2024') + 1]
    except (ValueError, IndexError):
        raise ValueError('Expected a HiForestpp2024 dataset path: ' + lfn)


def main():
    crab = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=crab / 'ppRefForestfiles',
                        help='Destination for lists and manifest (default: ppRefForestfiles)')
    args = parser.parse_args()
    source = crab / 'ppData2024ForestFiles.txt'
    mixing_pool = crab.parent / 'mixingFileList/ppRef2024_pool_pilot.txt'
    files = read_list(source)
    lfns = [logical_file_name(value) for value in files]
    if len(files) <= POOL_FILES or len(set(lfns)) != len(files):
        raise ValueError('Source needs at least 36 files and no duplicate LFNs')

    # Integer arithmetic makes the midpoint selection reproducible everywhere.
    indices = [(2 * index + 1) * len(files) // (2 * POOL_FILES)
               for index in range(POOL_FILES)]
    donors = [lfns[index] for index in indices]
    pool_lfns = set(donors)
    if len(pool_lfns) != POOL_FILES or not pool_lfns <= set(lfns):
        raise ValueError('Expected 35 unique donor files from the source list')
    pool_header = (
        '# ppRef 2024 mixing pilot: 35 reserved files from crab/ppData2024ForestFiles.txt\n'
        '# Source: supplied HardProbes-filtered forests, PPRefHardProbes0 through PPRefHardProbes4.\n'
        '# No independent ZeroBias or minimum-bias provenance is established for these donors.\n'
        '# TriggerSelection 8 selects signal events; the mixing reader disables donor trigger cuts.\n'
        '# Selection: evenly spaced source-list midpoints; see crab/ppRefForestfiles/manifest.json.\n'
    )
    pool_data = (pool_header + '\n'.join(REDIRECTOR + lfn for lfn in donors) + '\n').encode()
    signal = [value for value, lfn in zip(files, lfns) if lfn not in pool_lfns]

    outputs = {}
    tasks = []
    assigned_lfns = []
    for part, start in enumerate(range(0, len(signal), MAX_JOBS_PER_TASK), 1):
        entries = signal[start:start + MAX_JOBS_PER_TASK]
        # Retain the existing part/chunk selection interface; every part is one task.
        chunk = 1
        name = f'ppRefForestFiles_part{part:02d}_chunk{chunk:02d}.txt'
        data = ('\n'.join(entries) + '\n').encode()
        outputs[name] = data
        assigned_lfns.extend(logical_file_name(value) for value in entries)
        tasks.append({'part': part, 'chunk': chunk, 'file': name,
                      'files': len(entries), 'jobs': len(entries),
                      'sha256': sha256(data),
                      'request_name': f'eecAnalysis_ppRef_mixedPerp_part{part:02d}_chunk{chunk:02d}_Run3_v2'})
    signal_lfns = set(assigned_lfns)
    if (len(assigned_lfns) != len(signal_lfns) or signal_lfns & pool_lfns
            or signal_lfns | pool_lfns != set(lfns)):
        raise ValueError('Signal tasks must be disjoint and cover the source outside the donor pool')

    manifest = {
        'source': '../ppData2024ForestFiles.txt',
        'source_sha256': sha256(source.read_bytes()), 'source_files': len(files),
        'source_dataset_counts': dict(sorted(Counter(dataset_name(lfn) for lfn in lfns).items())),
        'mixing_pool': '../../mixingFileList/ppRef2024_pool_pilot.txt',
        'mixing_pool_sha256': sha256(pool_data), 'reserved_mixing_files': len(pool_lfns),
        'mixing_selection': {
            'strategy': '35 evenly spaced midpoints in the ordered source list',
            'index_formula': 'floor((2 * i + 1) * source_files / (2 * 35)), i = 0..34',
            'source_indices_zero_based': indices,
            'dataset_counts': dict(sorted(Counter(dataset_name(lfn) for lfn in donors).items())),
            'source_provenance': 'Supplied HardProbes-filtered forests; independent ZeroBias or minimum-bias provenance is not established',
            'trigger_policy': 'TriggerSelection 8 selects signal events; mixing reader disables donor trigger cuts',
            'event_content_validated': False,
        },
        'signal_files': len(signal), 'units_per_job': 1,
        'max_jobs_per_task': MAX_JOBS_PER_TASK, 'job_limit_source': LIMIT_SOURCE,
        'strategy': 'Consecutive source-order signal files after donor exclusion, packed into tasks of 10000 jobs with a smaller final task',
        'empty_file_policy': 'Retained; each empty signal forest affects only its own job',
        'tasks': tasks,
    }
    outputs['manifest.json'] = (json.dumps(manifest, indent=2) + '\n').encode()
    unexpected = {path.name for path in args.output_dir.glob('ppRefForestFiles_*.txt')} - set(outputs)
    if unexpected:
        raise ValueError('Archive previous lists before splitting: ' + ', '.join(sorted(unexpected)))
    targets = {args.output_dir / name: data for name, data in outputs.items()}
    targets[mixing_pool] = pool_data
    for path, data in targets.items():
        if path.exists() and path.read_bytes() != data:
            raise ValueError('Refusing to overwrite changed output: ' + str(path))
    for path, data in targets.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            with path.open('xb') as stream:
                stream.write(data)
    print(f'Validated {len(signal):,} signal files in {len(tasks)} tasks; '
          f'{len(pool_lfns)} mixing files excluded; one file per job.')
    for task in tasks:
        print(f"{task['file']}: {task['jobs']:,} jobs")


if __name__ == '__main__':
    main()
