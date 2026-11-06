#!/usr/bin/env python3
"""Split the ten original OO cohorts into tasks with one forest per job.

Preserve the original round-robin cohorts and reserved mixing pool. Each cohort
is subdivided in order at the CRAB default limit of 10,000 jobs per task.
Existing identical outputs are checked; changed outputs are never overwritten.
"""

import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit


PARTS = 10
MAX_JOBS_PER_TASK = 10000
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


def main():
    crab = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=crab / 'OOForestfiles',
                        help='Destination for lists and manifest (default: OOForestfiles)')
    args = parser.parse_args()
    source = crab / 'OOForestFiles.txt'
    mixing_pool = crab.parent / 'mixingFileList/OO2025_MB_pool_pilot.txt'
    files = read_list(source)
    lfns = [logical_file_name(value) for value in files]
    donors = read_list(mixing_pool)
    pool_lfns = {logical_file_name(value) for value in donors}
    if not files or len(set(lfns)) != len(files):
        raise ValueError('Source must be nonempty and contain no duplicate LFNs')
    if len(donors) != 35 or len(pool_lfns) != 35:
        raise ValueError('The OO pilot expects 35 unique mixing-pool files')
    if not pool_lfns <= set(lfns):
        raise ValueError('Reserved mixing-pool files are missing from the source list')
    signal = [value for value, lfn in zip(files, lfns) if lfn not in pool_lfns]
    outputs = {}
    tasks = []
    cohorts = []
    for part in range(1, PARTS + 1):
        cohort = signal[part - 1::PARTS]
        cohort_data = ('\n'.join(cohort) + '\n').encode()
        cohort_name = f'OOForestFiles_part{part:02d}.txt'
        archived = crab / 'OOForestfiles/previous_10parts' / cohort_name
        if archived.exists() and archived.read_bytes() != cohort_data:
            raise ValueError('Reconstructed cohort differs from archive: ' + str(archived))
        cohorts.append({'part': part, 'files': len(cohort),
                        'original_list': cohort_name, 'sha256': sha256(cohort_data)})
        for chunk, start in enumerate(range(0, len(cohort), MAX_JOBS_PER_TASK), 1):
            entries = cohort[start:start + MAX_JOBS_PER_TASK]
            name = f'OOForestFiles_part{part:02d}_chunk{chunk:02d}.txt'
            data = ('\n'.join(entries) + '\n').encode()
            outputs[name] = data
            tasks.append({'part': part, 'chunk': chunk, 'file': name,
                          'files': len(entries), 'jobs': len(entries),
                          'sha256': sha256(data),
                          'request_name': f'eecAnalysis_OO_mixedPerp_part{part:02d}_chunk{chunk:02d}_Run3_v2'})
    manifest = {
        'source': '../OOForestFiles.txt', 'source_sha256': sha256(source.read_bytes()),
        'source_files': len(files),
        'mixing_pool': '../../mixingFileList/OO2025_MB_pool_pilot.txt',
        'mixing_pool_sha256': sha256(mixing_pool.read_bytes()),
        'reserved_mixing_files': len(pool_lfns), 'signal_files': len(signal),
        'units_per_job': 1, 'max_jobs_per_task': MAX_JOBS_PER_TASK,
        'job_limit_source': LIMIT_SOURCE,
        'strategy': 'Original ten round-robin cohorts, then consecutive chunks of at most 10000 files',
        'empty_file_policy': 'Retained; each empty signal forest affects only its own job',
        'cohorts': cohorts, 'tasks': tasks,
    }
    outputs['manifest.json'] = (json.dumps(manifest, indent=2) + '\n').encode()
    unexpected = {path.name for path in args.output_dir.glob('OOForestFiles_*.txt')} - set(outputs)
    if unexpected:
        raise ValueError('Archive previous lists before splitting: ' + ', '.join(sorted(unexpected)))
    for name, data in outputs.items():
        path = args.output_dir / name
        if path.exists() and path.read_bytes() != data:
            raise ValueError('Refusing to overwrite changed output: ' + str(path))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, data in outputs.items():
        path = args.output_dir / name
        if not path.exists():
            with path.open('xb') as stream:
                stream.write(data)
    print(f'Validated {len(signal):,} signal files in {len(tasks)} tasks; '
          f'{len(pool_lfns)} mixing files excluded; one file per job.')
    for task in tasks:
        print(f"{task['file']}: {task['jobs']:,} jobs")


if __name__ == '__main__':
    main()
