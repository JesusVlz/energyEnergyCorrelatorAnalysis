#!/usr/bin/env python3
"""Checksum and retain closure ledgers, merge one frozen OO GenGen campaign."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

SUFFIXES = ('.root', '.partners.csv', '.pool.csv', '.jets.csv', '.contributions.csv')
PASS = 'PASS_FULL_FILE_SOFTWARE_AND_EXPOSURE_CHECKS'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def checked_file(path, expected):
    require(Path(path).stat().st_size == expected['size_bytes'], 'Size mismatch: ' + str(path))
    require(sha256(path) == expected['sha256'], 'SHA256 mismatch: ' + str(path))


def validate_record(record, entry, contract):
    """Bind the software PASS marker to the frozen production and complete publication."""
    require(record.get('status') == PASS, 'Unsuccessful run: ' + entry['stem'])
    for key in ('mode', 'data_type', 'archive_sha256', 'frozen_card_sha256',
                'pool_manifest_sha256', 'campaign_provenance_sha256', 'executable_sha256',
                'signal_headers_sha256', 'source_audit_sha256', 'partition_policy',
                'seed', 'pool_cap', 'holdout_policy'):
        require(record.get(key) == contract[key], 'Run contract mismatch: ' + key)
    for key in ('input_url', 'output_url'):
        require(record.get(key) == entry[key], 'Run URL mismatch: ' + key)
    require(record['campaign_provenance'] == contract['campaign_provenance'], 'Campaign provenance mismatch')
    source = record['source']
    for key in ('source_entries', 'root_uuid', 'size_bytes'):
        require(source[key] == entry['source'][key], 'Source identity mismatch: ' + key)
    require(record['source_identity']['original_path'] == entry['source']['path'], 'Source path mismatch')
    require(record['source_identity']['root_uuid'] == source['root_uuid'], 'Source identity UUID mismatch')
    require(record['source_identity']['full_event_entries'] == source['source_entries'], 'Truncated source')
    require(record['root_validation']['all_events'] == source['source_entries'], 'ROOT event exposure mismatch')
    require(record['ledger_validation']['source_entries'] == source['source_entries'], 'Ledger source exposure mismatch')
    require(record['root_validation']['inclusive_jet_entries'] == record['ledger_validation']['matched_jets'], 'Jet exposure mismatch')
    require(record['ledger_validation']['pool_entries'] == contract['pool_cap'], 'Incomplete donor scan')
    require(record.get('commands') and all(command['exit_code'] == 0 for command in record['commands']), 'Failed worker command')
    names = {entry['stem'] + suffix for suffix in SUFFIXES}
    require(set(record['artifacts']) == names, 'Incomplete or unexpected artifact set')
    urls = {entry['output_url'][:-5] + suffix for suffix in SUFFIXES}
    require(set(record['published_artifacts']) == urls, 'Incomplete publication')
    for name, item in record['artifacts'].items():
        require(item['size_bytes'] > 0 and len(item['sha256']) == 64
                and all(char in '0123456789abcdef' for char in item['sha256']), 'Invalid artifact checksum: ' + name)


def validate_manifest(manifest):
    require(manifest.get('schema_version') == 1 and manifest.get('kind') == 'OOMC_GenGen_closure_merge', 'Unsupported merge manifest')
    require(manifest['contract']['mode'] == 'GenGen' and manifest['contract']['data_type'] == 10, 'Unsupported merge mode')
    require(len(manifest['entries']) == manifest['input_count'] > 0, 'Incomplete merge inputs')
    for field in ('stem', 'input_url', 'output_url'):
        values = [entry[field] for entry in manifest['entries']]
        require(len(set(values)) == len(values), 'Duplicate merge input: ' + field)
    for entry in manifest['entries']:
        require(Path(entry['stem']).name == entry['stem'] and entry['stem'] not in ('.', '..'), 'Unsafe output stem')
        require(entry['output_url'].startswith(manifest['source_output_url'] + '/'), 'Mixed production outputs')
        validate_record(entry['run_record'], entry, manifest['contract'])
    for key in ('all_events', 'inclusive_jet_entries', 'eec_entries'):
        total = sum(entry['run_record']['root_validation'][key] for entry in manifest['entries'])
        require(manifest['expected_root_totals'][key] == total, 'Manifest exposure total mismatch: ' + key)


def run(manifest_path, work_dir, publish=True):
    from merge_root import merge_root
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    validate_manifest(manifest)
    work = Path(work_dir).resolve()
    require(not work.exists(), 'Work directory already exists')
    work.mkdir(parents=True)
    commands = []

    def command(argv):
        result = subprocess.run([str(value) for value in argv], check=False)
        commands.append({'argv': [str(value) for value in argv], 'exit_code': result.returncode})
        require(result.returncode == 0, 'Worker command failed: ' + str(argv))

    inputs = []
    bundle = work / 'merged_OOMC_GenGen.ledgers.tar.gz'
    # Keep original per-job ledgers: repeated pools retain their job-specific reuse counts.
    with tarfile.open(bundle, 'w:gz', compresslevel=1) as archive:
        archive.add(manifest_path, arcname='merge_manifest.json', recursive=False)
        for index, entry in enumerate(manifest['entries'], 1):
            stem = entry['stem']
            local_record = work / (stem + '.run.json')
            command(['xrdcp', '--nopbar', '--retry', '3', '--cksum', 'adler32',
                     entry['output_url'][:-5] + '.run.json', local_record])
            checked_file(local_record, entry['run_metadata'])
            record = json.loads(local_record.read_text())
            validate_record(record, entry, manifest['contract'])
            archive.add(local_record, arcname=local_record.name, recursive=False)
            local_record.unlink()
            for suffix in SUFFIXES:
                name = stem + suffix
                local = work / name
                command(['xrdcp', '--nopbar', '--retry', '3', '--cksum', 'adler32',
                         entry['output_url'][:-5] + suffix, local])
                checked_file(local, record['artifacts'][name])
                if suffix == '.root':
                    inputs.append(local)
                else:
                    archive.add(local, arcname=name, recursive=False)
                    local.unlink()
            print('MERGE_INPUT_VERIFIED ' + str(index) + '/' + str(len(manifest['entries'])), flush=True)
    output = work / 'merged_OOMC_GenGen.root'
    validation = merge_root(inputs, output, manifest['expected_root_totals'])
    report = {'schema_version': 1, 'status': 'PASS_MERGE_SOFTWARE_AND_EXPOSURE_CHECKS',
              'campaign': manifest['tag'], 'input_count': len(inputs),
              'manifest_sha256': sha256(manifest_path), 'root_validation': validation,
              'finished_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'limitations': manifest['limitations'], 'artifacts': {}}
    artifacts = [output, bundle]
    for path in artifacts:
        report['artifacts'][path.name] = {'sha256': sha256(path), 'size_bytes': path.stat().st_size}
    report_path = work / 'merged_OOMC_GenGen.run.json'
    if publish:
        for path in artifacts:
            # Existing results require review; never silently replace them.
            command(['xrdcp', '--nopbar', '--retry', '3', '--cksum', 'adler32', '--posc',
                     path, manifest['destination_url'] + '/' + path.name])
    report['commands'] = commands
    report['published'] = publish
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    if publish:
        command(['xrdcp', '--nopbar', '--retry', '3', '--cksum', 'adler32', '--posc',
                 report_path, manifest['destination_url'] + '/' + report_path.name])
    print('MERGE_JOB_SUCCESS ' + json.dumps({'campaign': manifest['tag'], 'inputs': len(inputs), 'published': publish}), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--no-publish', action='store_true')
    args = parser.parse_args()
    run(args.manifest, args.work_dir, not args.no_publish)


if __name__ == '__main__':
    main()
