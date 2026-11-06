#!/usr/bin/env python3
"""Run the established pair and response helpers on one staged forest.

All paths are worker-local. This driver neither transfers files nor submits jobs.
Completion markers are emitted only after both helpers and their shared physical
event ledger pass the full-file coverage and provenance checks.
"""
import argparse
from collections import deque
import csv
import gzip
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import time
from urllib.parse import urlsplit

from frozen_input import stage_frozen_input

IDENTITY = struct.Struct('<QQQB')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def load(path):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as stream:
        return json.load(stream)


def save(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def checked_artifact(directory, record):
    name = Path(record['name'])
    if name.is_absolute() or '..' in name.parts:
        raise ValueError('Unsafe pair artifact name')
    path = directory/name
    if path.stat().st_size != record['bytes'] or digest(path) != record['sha256']:
        raise ValueError('Pair artifact checksum mismatch: '+str(name))
    return path


def compare_ledgers(pair_ledger, response_ledger, expected, selected_counts):
    counts = {'raw': 0, 'train': 0, 'holdout': 0}
    with gzip.open(pair_ledger, 'rb') as pairs, gzip.open(response_ledger, 'rt') as response:
        for index, row in enumerate(csv.DictReader(response)):
            record = pairs.read(IDENTITY.size)
            if len(record) != IDENTITY.size:
                raise ValueError('Pair/response raw identity coverage differs')
            run, lumi, event, marker = IDENTITY.unpack(record)
            physical = tuple(int(row[key]) for key in ('run', 'lumi', 'event'))
            if int(row['source_entry']) != index or (run, lumi, event) != physical:
                raise ValueError('Pair/response physical identity mismatch at source entry '+str(index))
            fold = int(row['fold'])
            if fold not in (0, 1) or marker not in (0, 1, 2) or (marker and marker != fold+1):
                raise ValueError('Pair/response training-validation split mismatch')
            counts['raw'] += 1
            if marker:
                counts['train' if marker == 1 else 'holdout'] += 1
        if pairs.read(1) or counts['raw'] != expected:
            raise ValueError('Pair/response ledger does not cover the full frozen input')
    if any(counts[key] != selected_counts[key] for key in ('raw', 'train', 'holdout')):
        raise ValueError('Pair selected ledger/receipt counts disagree')
    return {'status': 'PASS', 'matched_raw_physical_ids': counts['raw'],
            'selected_train': counts['train'], 'selected_holdout': counts['holdout']}


def write_marker(directory, provenance):
    files = sorted(path for path in directory.rglob('*') if path.is_file() and path.name != 'stage_complete.json')
    marker = {'status': 'COMPLETE', 'files': {
        str(path.relative_to(directory)): {'bytes': path.stat().st_size, 'sha256': digest(path)}
        for path in files}, 'frozen_task': provenance['frozen_task'],
        'campaign_manifest_sha256': provenance['campaign_manifest_sha256']}
    save(directory/'stage_complete.json', marker)


def failure_tails(output, error):
    print('PAIR_RESPONSE_WORKER_FAILED: '+repr(error), file=sys.stderr, flush=True)
    for path in sorted(output.rglob('*.log')):
        try:
            with path.open(errors='replace') as stream:
                tail = ''.join(deque(stream, maxlen=80))
            print('Native log tail: '+str(path)+'\n'+tail, file=sys.stderr, flush=True)
        except OSError as problem:
            print('Cannot read native log: '+str(problem), file=sys.stderr, flush=True)


def run(args):
    analysis = args.analysis_dir.resolve(strict=True)
    manifest_path = args.manifest.resolve(strict=True)
    source = args.input_local.resolve()
    output = args.output_dir.resolve()
    scratch = args.scratch.resolve()
    manifest = load(manifest_path)
    if not 0 <= args.task_index < len(manifest['tasks']):
        raise ValueError('Task index outside frozen manifest')
    task = manifest['tasks'][args.task_index]
    expected = task['entry_limit']
    identity = task['input_identity']
    if expected <= 0 or identity.get('entries', expected) != expected:
        raise ValueError('Manifest task is not a complete nonempty input file')
    source_metadata = None
    input_url = getattr(args, 'input_url', None)
    if input_url:
        actual_path = '/' + urlsplit(input_url).path.lstrip('/')
        expected_path = str(task['input']).replace('/eos/home-j/', '/eos/user/j/')
        if actual_path != expected_path:
            raise ValueError('Input URL differs from frozen task source')
        if source.exists():
            raise ValueError('Input staging destination already exists')
        source.parent.mkdir(parents=True, exist_ok=True)
        source_metadata = stage_frozen_input(input_url, source, identity,
            lambda: subprocess.run(['xrdcp', '--retry', '3', '--cksum', 'adler32',
                                    input_url, str(source)], check=True))
    elif not identity.get('sha256'):
        raise ValueError('Already staged input requires a frozen SHA256; use --input-url to verify source size/mtime while staging')
    if source.stat().st_size != identity['size_bytes']:
        raise ValueError('Staged input size differs from frozen input')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Combined output directory must be new or empty')
    output.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    pairs, response = output/'pairs', output/'response'
    pair_helper = analysis/'validation/run3Pairs/batch_worker.py'
    response_helper = analysis/'validation/run3Unfolding/batch/response_worker.py'
    fingerprint_paths = [analysis/'eecAnalysis', pair_helper, response_helper,
        analysis/'validation/run3Pairs/pair_counts.py', analysis/'validation/run3Pairs/batch_pairs.py',
        analysis/'validation/run3Pairs/export_forest.py',
        response_helper.parent/'cards'/(task['system']+'_R04_response.input')]
    fingerprints = {str(path.relative_to(analysis)): digest(path) for path in fingerprint_paths}
    source_digest = digest(source)
    if 'sha256' in identity and identity['sha256'] != source_digest:
        raise ValueError('Staged input checksum differs from frozen input')
    provenance = {'frozen_task': task, 'campaign_manifest_sha256': digest(manifest_path),
        'source_input_sha256': source_digest, 'source_metadata': source_metadata, 'source_fingerprints': fingerprints,
        'driver_sha256': digest(__file__), 'task_index': args.task_index,
        'scope': 'Pair training shards and independent native response matrices; no unfolding'}
    started = time.monotonic()
    try:
        subprocess.run([sys.executable, str(pair_helper), '--manifest', str(manifest_path),
            '--task-index', str(args.task_index), '--input-local', str(source),
            '--output-dir', str(pairs), '--scratch', str(scratch)], cwd=analysis, check=True)
        prefix = '{}_{:03d}'.format(task['system'], task['file_index'])
        pair_receipt = load(pairs/(prefix+'.receipt.json'))
        if pair_receipt['task'] != task or pair_receipt['event_counts']['raw'] != expected:
            raise ValueError('Pair helper did not cover the frozen full-file task')
        artifacts = {key: checked_artifact(pairs, record) for key, record in pair_receipt['files'].items()}
        counts = load(artifacts['counts'])
        metadata = counts['native_metadata']
        if counts['task'] != task or counts['event_counts'] != pair_receipt['event_counts']:
            raise ValueError('Pair counts/receipt task mismatch')
        if counts['raw_identity_scope'] != 'ALL_SCANNED_EVENTS':
            raise ValueError('Pair ledger is not a full raw event ledger')
        if (metadata['raw_input_sha256'] != source_digest or
            metadata['frozen_task_manifest_sha256'] != provenance['campaign_manifest_sha256'] or
            metadata['native_binary_sha256'] != fingerprints['eecAnalysis'] or
            counts['physics_core_sha256'] != fingerprints['validation/run3Pairs/pair_counts.py']):
            raise ValueError('Pair source/native/core fingerprint mismatch')
        subprocess.run([sys.executable, str(response_helper), str(source), str(response),
            '--system', task['system'], '--namespace', task['sample_namespace'],
            '--source-input', task['input'], '--analysis-repo', str(analysis),
            '--scratch-parent', str(scratch)], cwd=analysis, check=True)
        response_receipt = load(response/'receipt.json')
        skim = response_receipt['skim']
        if (response_receipt['status'] != 'COMPLETE_ENGINEERING_ONLY' or
            skim['input_entries'] != expected or skim['processed_entries'] != expected or
            sum(skim['fold_entries']) != expected or skim['sample_namespace'] != task['sample_namespace']):
            raise ValueError('Response does not cover every frozen input entry')
        if (response_receipt['native_sha256'] != fingerprints['eecAnalysis'] or
            response_receipt['pair_counts_sha256'] != fingerprints['validation/run3Pairs/pair_counts.py'] or
            response_receipt['helper_sha256'] != fingerprints[str(response_helper.relative_to(analysis))] or
            response_receipt['card_sha256'] != fingerprints[str(fingerprint_paths[-1].relative_to(analysis))]):
            raise ValueError('Response source/native/card fingerprint mismatch')
        for label in ('train', 'validation'):
            report = response_receipt['folds'][label]
            path = response/(label+'_raw.root')
            if path.stat().st_size != report['bytes'] or digest(path) != report['sha256']:
                raise ValueError('Response output checksum mismatch: '+label)
        if digest(response/'event_identity.csv.gz') != skim['ledger_sha256']:
            raise ValueError('Response identity ledger checksum mismatch')
        provenance['shared_identity_audit'] = compare_ledgers(artifacts['identities'], response/'event_identity.csv.gz', expected, pair_receipt['event_counts'])
        response_receipt.update(frozen_task=task, campaign_manifest_sha256=provenance['campaign_manifest_sha256'])
        save(response/'receipt.json', response_receipt)
        provenance.update(status='COMPLETE', elapsed_seconds=time.monotonic()-started)
        for directory in (pairs, response):
            save(directory/'combined_worker_provenance.json', provenance)
        for directory in (pairs, response):
            write_marker(directory, provenance)
    except Exception as error:
        provenance.update(status='FAILED', error=repr(error), elapsed_seconds=time.monotonic()-started)
        save(output/'failure.json', provenance)
        failure_tails(output, error)
        raise
    print(json.dumps({'status': 'PAIR_RESPONSE_LOCAL_COMPLETE', 'output': str(output),
        'raw_entries': expected, 'task_index': args.task_index}), flush=True)
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('analysis-dir', 'manifest', 'input-local', 'output-dir', 'scratch'):
        parser.add_argument('--'+flag, type=Path, required=True)
    parser.add_argument('--input-url', help='Stage the frozen source after checking remote size and modification time')
    parser.add_argument('--task-index', type=int, required=True)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
