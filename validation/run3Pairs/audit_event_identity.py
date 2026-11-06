#!/usr/bin/env python3
"""Refresh the existing raw-input event-identity audit for completed pair exports."""
import argparse
import hashlib
import json
import shutil
import time
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(manifest_path):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    manifest = json.loads(manifest_path.read_text())
    systems, identities, physical = {}, {}, set()
    files = []
    for task in manifest['tasks']:
        system = task['system']
        limit = int(task.get('entry_limit', manifest.get('events_per_file')))
        metadata_path = Path(task['output'] + '.metadata.json')
        metadata = json.loads(metadata_path.read_text())
        namespace = metadata['sample_namespace']
        if metadata['entry_limit_per_file'] != limit:
            raise ValueError('Metadata entry limit differs from finalized manifest')
        source = Path(task['input'])
        arrays = ROOT.RDataFrame('hiEvtAnalyzer/HiTree', str(source)).Range(limit).AsNumpy(['run', 'lumi', 'evt'])
        count = len(arrays['evt'])
        if count != limit:
            raise ValueError(f'{source}: expected {limit} scanned entries, found {count}')
        ledger_path = Path(metadata['native_outputs'])/'reco.jsonl.events.jsonl'
        seen = identities.setdefault(system, set())
        summary = systems.setdefault(system, {'files': 0, 'scanned_rows': 0, 'unique_physical_ids': 0,
                                              'duplicate_rows': 0, 'examples': [], 'raw_native_ledger_mismatches': 0})
        summary['files'] += 1
        duplicates, mismatches = 0, 0
        with ledger_path.open() as ledger:
            for index, (run, lumi, evt) in enumerate(zip(arrays['run'], arrays['lumi'], arrays['evt'])):
                identity = (int(run), int(lumi), int(evt))
                row = ledger.readline()
                if not row:
                    raise ValueError(f'{ledger_path}: ledger ends before raw source entry {index}')
                native = json.loads(row)
                if native['source_entry'] != index or tuple(native['id']) != identity:
                    mismatches += 1
                if identity in seen:
                    duplicates += 1
                    if len(summary['examples']) < 10:
                        summary['examples'].append({'source': str(source), 'entry': index, 'id': list(identity)})
                seen.add(identity)
                physical.add((namespace, *identity))
            if ledger.readline():
                raise ValueError(f'{ledger_path}: ledger has more than {limit} raw entries')
        summary['scanned_rows'] += count
        summary['duplicate_rows'] += duplicates
        summary['raw_native_ledger_mismatches'] += mismatches
        summary['unique_physical_ids'] = len(seen)
        files.append({'system': system, 'sample_namespace': namespace, 'input': str(source),
                      'raw_file_size_bytes': source.stat().st_size, 'entry_range': [0, count],
                      'scanned_rows': count, 'duplicates_within_system': duplicates,
                      'raw_native_ledger_mismatches': mismatches, 'native_ledger': str(ledger_path),
                      'native_ledger_sha256': digest(ledger_path)})
        print(json.dumps({'system': system, 'file': source.name, 'scanned_rows': count,
                          'duplicate_rows': duplicates, 'ledger_mismatches': mismatches}), flush=True)
    total = sum(row['scanned_rows'] for row in systems.values())
    passed = all(not row['duplicate_rows'] and not row['raw_native_ledger_mismatches'] for row in systems.values())
    return {**systems, 'status': 'PASS' if passed and len(physical) == total else 'FAIL',
            'global': {'files': len(files), 'scanned_rows': total, 'unique_namespaced_physical_ids': len(physical),
                       'duplicate_namespaced_rows': total-len(physical)},
            'method': 'Read exact integer run/lumi/evt from hiEvtAnalyzer/HiTree in original ROOT sources for the manifest entry range; verify every row against final native scanned-event ledger. Namespace distinguishes independent collision samples.',
            'input_manifest': str(manifest_path.resolve()), 'input_manifest_sha256': digest(manifest_path),
            'tool_sha256': digest(Path(__file__)), 'files': files}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--backup', type=Path, required=True, help='preserve previous audit here without overwriting it')
    args = parser.parse_args()
    started = time.monotonic()
    if args.output.exists() and not args.backup.exists():
        shutil.copy2(args.output, args.backup)
    result = audit(args.manifest)
    result['previous_audit'] = ({'path': str(args.backup.resolve()), 'sha256': digest(args.backup)}
                                if args.backup.exists() else None)
    result['elapsed_seconds'] = time.monotonic()-started
    temporary = args.output.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    temporary.replace(args.output)
    print(json.dumps({'status': result['status'], 'global': result['global'],
                      'systems': {system: result[system] for system in ('OO', 'pp')}, 'elapsed_seconds': result['elapsed_seconds']}))
    if result['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
