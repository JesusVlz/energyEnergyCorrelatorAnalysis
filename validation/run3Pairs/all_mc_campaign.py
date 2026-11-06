#!/usr/bin/env python3
"""Prepare, run, or validate a frozen full-file Run-3 R4 MC export campaign.

prepare reads ROOT headers only. run is the only command that launches exports.
A new campaign never reuses exports from another campaign/source snapshot.
"""
import argparse
import concurrent.futures
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

from run_sample import REPO, SAMPLES, digest
from export_forest import PAYLOADS

VERSION = 'run3-an751-all-mc-v1'
CARDS = {'OO': 'cardEECOOMCInclusiveR04.input', 'pp': 'cardEECpp5v36TeVMCR04.input'}


def write_json(path, data):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def identity(path, hashed=False):
    path = path.resolve(strict=True)
    stat = path.stat()
    result = dict(path=str(path), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if hashed:
        result['sha256'] = digest(path)
    return result


def verify_identity(expected, hashed=False):
    actual = identity(Path(expected['path']), hashed)
    for key in actual:
        if actual[key] != expected[key]:
            raise ValueError(f"Changed frozen {key}: {expected['path']}")


def root_entries(path):
    import ROOT
    source = ROOT.TFile.Open(str(path), 'READ')
    if not source:
        raise OSError(f'Cannot open ROOT input: {path}')
    try:
        if source.IsZombie():
            raise OSError(f'Invalid ROOT input: {path}')
        tree = source.Get('hiEvtAnalyzer/HiTree')
        if not tree or not tree.InheritsFrom('TTree'):
            raise ValueError(f'Missing hiEvtAnalyzer/HiTree: {path}')
        return int(tree.GetEntries())
    finally:
        source.Close()


def discover(directory, prefix):
    """Use actual numeric matches: a missing index does not truncate discovery."""
    found = []
    for path in Path(directory).glob(prefix + '*.root'):
        match = re.fullmatch(re.escape(prefix) + r'(\d+)\.root', path.name)
        if not match:
            raise ValueError(f'Unexpected matching filename: {path}')
        found.append((int(match[1]), path.resolve(strict=True)))
    if not found or len({number for number, _ in found}) != len(found):
        raise ValueError(f'Empty or ambiguous input discovery: {directory}/{prefix}*.root')
    return sorted(found)


def production_sources():
    """Freeze native source/binary, production scripts, cards and active payloads."""
    paths = {REPO/'eecAnalysis', REPO/'eecAnalysis.cxx', REPO/'Makefile'}
    paths.update((REPO/'src').glob('*.h'))
    paths.update((REPO/'src').glob('*.cxx'))
    paths.update(REPO/'validation/run3Pairs'/name for name in (
        'all_mc_campaign.py', 'run_sample.py', 'export_forest.py', 'pair_counts.py',
        'plot_closure.py', 'summarize_results.py', 'audit_event_identity.py'))
    paths.update(PAYLOADS.values())
    paths.update(REPO/'jetEnergyCorrections'/name for name in (
        'L2Relative_AK4PF_OO_V1.txt', 'Prompt24HIpp_V1_MC_L2Relative_AK4PF.txt',
        'CorrectionNotAppliedPF.txt'))
    for card_name in CARDS.values():
        card = REPO/card_name
        paths.add(card)
        for line in card.read_text().splitlines():
            fields = line.split('#', 1)[0].split()
            if fields and fields[0] in ('JetVetoMapFile', 'Run3TrackPairCorrectionFile'):
                paths.add(REPO/fields[1])
    return sorted(path.resolve(strict=True) for path in paths)


def prepare(output, samples=None, sources=None, entry_reader=root_entries, repo=REPO):
    samples = SAMPLES if samples is None else samples
    sources = production_sources() if sources is None else sources
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f'Prepare requires a new or empty output directory: {output}')
    snapshot = [identity(path, hashed=True) for path in sources]
    tasks, discovery = [], {}
    for system, (directory, prefix, namespace) in samples.items():
        found = discover(directory, prefix)
        indices = [number for number, _ in found]
        discovery[system] = dict(directory=str(Path(directory).resolve()), pattern=prefix+'*.root',
                                 files=len(found), indices=indices,
                                 missing_indices=[i for i in range(1, max(indices)+1) if i not in set(indices)])
        for number, path in found:
            frozen = identity(path)
            entries = entry_reader(path)
            verify_identity(frozen)
            if entries < 1:
                raise ValueError(f'Empty raw event tree; cannot use positive full-file bound: {path}')
            target = output/f'{system}_{number:03d}.jsonl'
            command = [str(Path(sys.executable).resolve()), str(repo/'validation/run3Pairs/export_forest.py'),
                       str(path), str(target), '--system', system, '--sample-namespace', namespace,
                       '--max-events', str(entries)]
            tasks.append(dict(system=system, file_index=number, sample_namespace=namespace,
                              input=str(path), output=str(target), entry_limit=entries,
                              input_identity={**frozen, 'entries': entries}, command=command))
            if len(tasks) % 50 == 0:
                print(json.dumps({'prepared_headers': len(tasks), 'last_input': str(path)}), flush=True)
    for source in snapshot:
        verify_identity(source, hashed=True)
    counts = {system: {'files': sum(t['system'] == system for t in tasks),
                      'raw_entries': sum(t['entry_limit'] for t in tasks if t['system'] == system)}
              for system in samples}
    bounds = {task['entry_limit'] for task in tasks}
    manifest = dict(procedure=VERSION, created_utc=datetime.now(timezone.utc).isoformat(),
                    repository=str(repo.resolve()), python_executable=str(Path(sys.executable).resolve()),
                    tasks=tasks, discovery=discovery, expected=counts,
                    events_per_file=next(iter(bounds)) if len(bounds) == 1 else None,
                    entry_policy='Each task scans every raw hiEvtAnalyzer/HiTree entry frozen at prepare.',
                    training_validation='unchanged namespaced physical-event FNV1a64 modulo2',
                    source_identity=snapshot,
                    source_sha256={str(Path(p['path']).relative_to(repo.resolve())): p['sha256'] for p in snapshot})
    output.mkdir(parents=True, exist_ok=True)
    write_json(output/'sample_manifest.json', manifest)
    write_json(output/'frozen_commands.json', [t['command'] for t in tasks])
    (output/'frozen_commands.txt').write_text('\n'.join(shlex.join(t['command']) for t in tasks)+'\n')
    for system in ('all', *samples):
        name = 'inputs.txt' if system == 'all' else f'{system}_inputs.txt'
        (output/name).write_text('\n'.join(t['output'] for t in tasks if system == 'all' or t['system'] == system)+'\n')
    print(json.dumps({'status': 'PREPARED_ONLY', 'output': str(output), 'expected': counts,
                      'manifest_sha256': digest(output/'sample_manifest.json')}), flush=True)
    return manifest


def load_manifest(output):
    path = output.resolve()/'sample_manifest.json'
    manifest = json.loads(path.read_text())
    if manifest['procedure'] != VERSION:
        raise ValueError('Not a frozen all-MC campaign manifest')
    return manifest, digest(path)


def verify_frozen(manifest, tasks, headers=False, full_sources=True):
    for source in manifest['source_identity']:
        verify_identity(source, hashed=full_sources)
    for task in tasks:
        verify_identity(task['input_identity'])
        if headers and root_entries(Path(task['input'])) != task['entry_limit']:
            raise ValueError('Changed raw tree entry count: '+task['input'])


def artifact_scan(task, manifest):
    """Validate metadata and count the native scanned-entry ledger, not selected rows."""
    target = Path(task['output'])
    metadata_path = Path(str(target)+'.metadata.json')
    metadata = json.loads(metadata_path.read_text())
    for key, expected in [('system', task['system']), ('input', task['input']),
                          ('sample_namespace', task['sample_namespace']),
                          ('entry_limit_per_file', task['entry_limit']),
                          ('input_size_bytes', task['input_identity']['size_bytes'])]:
        if metadata.get(key) != expected:
            raise ValueError(f'Export metadata mismatch ({key}): {target}')
    hashes = manifest['source_sha256']
    for key, source in [('tool_sha256', 'validation/run3Pairs/export_forest.py'),
                        ('analyzer_sha256', 'src/EECAnalyzer.cxx')]:
        if metadata.get(key) != hashes[source]:
            raise ValueError(f'Export source mismatch ({key}): {target}')
    payload = str(Path(metadata['payload']).resolve())
    frozen_payload = next((p for p in manifest['source_identity'] if p['path'] == payload), None)
    if frozen_payload is None or metadata['payload_sha256'] != frozen_payload['sha256']:
        raise ValueError(f'Export calibration is outside frozen snapshot: {target}')
    output_hash = digest(target)
    if output_hash != metadata['output_sha256']:
        raise ValueError(f'Export SHA256 mismatch: {target}')
    work = target.parent/(target.stem+'.native')
    if Path(metadata['native_outputs']).resolve() != work:
        raise ValueError(f'Unexpected native artifact directory: {target}')
    artifacts = {str(target): output_hash, str(metadata_path): digest(metadata_path)}
    ledger_hashes, scanned = [], []
    for population in ('reco', 'gen'):
        ledger = work/(population+'.jsonl.events.jsonl')
        count, checksum = 0, hashlib.sha256()
        with ledger.open('rb') as stream:
            for line in stream:
                checksum.update(line)
                row = json.loads(line)
                if row['source_entry'] != count:
                    raise ValueError(f'Nonsequential native scanned-event ledger: {ledger}')
                count += 1
        scanned.append(count)
        ledger_hashes.append(checksum.hexdigest())
        artifacts[str(ledger)] = ledger_hashes[-1]
        for suffix in ('.root', '.input', '.jsonl', '.log'):
            native = work/(population+suffix)
            artifacts[str(native)] = digest(native)
    if scanned != [task['entry_limit']]*2 or ledger_hashes[0] != ledger_hashes[1]:
        raise ValueError(f'Native ledgers do not match full frozen input entries: {target}')
    return dict(metadata=metadata, artifacts_sha256=artifacts, scanned_entries=scanned[0])


def receipt_path(task):
    return Path(task['output']).with_suffix('.run.json')


def validate_success(task, manifest, campaign_hash):
    saved = json.loads(receipt_path(task).read_text())
    if saved.get('campaign_sha256') != campaign_hash or saved.get('command') != task['command']:
        raise ValueError('Receipt is from another campaign/command: '+task['output'])
    if saved.get('exit_code') != 0:
        return saved
    actual = artifact_scan(task, manifest)
    if actual != {key: saved.get(key) for key in actual}:
        raise ValueError('Successful export or metadata changed: '+task['output'])
    verify_identity(task['input_identity'])
    if digest(Path(task['input'])) != saved.get('input_sha256'):
        raise ValueError('Successful export input changed: '+task['input'])
    return saved


def has_artifacts(task):
    target = Path(task['output'])
    return any((target.exists(), Path(str(target)+'.metadata.json').exists(),
                target.with_suffix('.export.log').exists(),
                (target.parent/(target.stem+'.native')).exists()))


def archive_failure(task):
    target = Path(task['output'])
    archive = target.parent/'failures'/target.stem/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive.mkdir(parents=True)
    for path in (receipt_path(task), target, Path(str(target)+'.metadata.json'),
                 target.with_suffix('.export.log'), target.parent/(target.stem+'.native')):
        if path.exists():
            shutil.move(str(path), str(archive/path.name))


def run_task(task, manifest, campaign_hash):
    verify_frozen(manifest, [task], full_sources=False)
    if receipt_path(task).exists():
        archive_failure(task)
    started = time.monotonic()
    result = dict(task, campaign_sha256=campaign_hash, source_sha256=manifest['source_sha256'], exit_code=1)
    try:
        with Path(task['output']).with_suffix('.export.log').open('w') as log:
            result['exit_code'] = subprocess.run(task['command'], cwd=manifest['repository'],
                                                  stdout=log, stderr=subprocess.STDOUT).returncode
        if result['exit_code'] == 0:
            verify_frozen(manifest, [task], full_sources=False)
            result.update(artifact_scan(task, manifest))
            result['input_sha256'] = digest(Path(task['input']))
            verify_identity(task['input_identity'])
    except Exception as error:
        result.update(exit_code=1, error=f'{type(error).__name__}: {error}')
    result['elapsed_seconds'] = time.monotonic()-started
    write_json(receipt_path(task), result)
    print(json.dumps({key: result[key] for key in ('system', 'input', 'exit_code', 'elapsed_seconds')}), flush=True)
    return result


def summarize(manifest, results):
    systems = {}
    for system in manifest['expected']:
        successes = [r for r in results if r['system'] == system and r['exit_code'] == 0]
        counts = {}
        for result in successes:
            for key, value in result['metadata']['counts'].items():
                counts[key] = counts.get(key, 0)+value
        systems[system] = dict(expected=manifest['expected'][system], completed_files=len(successes),
                               scanned_entries=sum(r['scanned_entries'] for r in successes), counts=counts)
    failures = [r['input'] for r in results if r['exit_code'] != 0]
    complete = all(s['completed_files'] == s['expected']['files'] and
                   s['scanned_entries'] == s['expected']['raw_entries'] for s in systems.values())
    return dict(status='COMPLETE' if complete and not failures else 'INCOMPLETE', systems=systems, failed=failures)


@contextmanager
def campaign_lock(output):
    with (output/'.campaign.lock').open('a') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError('Another campaign run/status holds the output lock') from error
        yield


def execute(output, workers=1, resume=False, system='all'):
    manifest, campaign_hash = load_manifest(output)
    selected = [t for t in manifest['tasks'] if system == 'all' or t['system'] == system]
    verify_frozen(manifest, selected, headers=True)
    pending = []
    # Complete all preflight checks before launching any subprocess.
    for task in selected:
        if receipt_path(task).exists():
            if not resume:
                raise ValueError('Existing receipt requires --resume: '+str(receipt_path(task)))
            saved = validate_success(task, manifest, campaign_hash)
            if saved['exit_code'] == 0:
                continue
        elif has_artifacts(task):
            raise ValueError('Unreceipted partial artifacts require inspection: '+task['output'])
        pending.append(task)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_task, task, manifest, campaign_hash): task for task in pending}
        for future in concurrent.futures.as_completed(futures):
            if future.cancelled():
                continue
            try:
                failed = future.result()['exit_code'] != 0
            except Exception:
                for queued in futures:
                    queued.cancel()
                raise
            if failed:
                for queued in futures:
                    queued.cancel()
    verify_frozen(manifest, selected)
    return status(output)


def status(output):
    manifest, campaign_hash = load_manifest(output)
    verify_frozen(manifest, manifest['tasks'])
    results, pending = [], []
    for task in manifest['tasks']:
        if receipt_path(task).exists():
            results.append(validate_success(task, manifest, campaign_hash))
        else:
            if has_artifacts(task):
                raise ValueError('Unreceipted partial artifacts: '+task['output'])
            pending.append(task['input'])
    summary = {**summarize(manifest, results), 'campaign_sha256': campaign_hash,
               'pending_files': len(pending), 'source_sha256': manifest['source_sha256']}
    write_json(output/'sample_summary.json', summary)
    print(json.dumps({k: summary[k] for k in ('status', 'systems', 'pending_files', 'failed')}), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    for action in ('prepare', 'run', 'status'):
        command = sub.add_parser(action)
        command.add_argument('output', type=Path)
        if action == 'run':
            command.add_argument('--workers', type=int, choices=range(1, 5), default=1)
            command.add_argument('--resume', action='store_true')
            command.add_argument('--system', choices=('OO', 'pp', 'all'), default='all')
    args = parser.parse_args()
    output = args.output.resolve()
    if args.action == 'prepare':
        prepare(output)
    else:
        with campaign_lock(output):
            summary = execute(output, args.workers, args.resume, args.system) if args.action == 'run' else status(output)
        if summary['failed']:
            raise SystemExit(1)
        if summary['status'] != 'COMPLETE' and (args.action == 'status' or args.system == 'all'):
            raise SystemExit(2)


if __name__ == '__main__':
    main()
