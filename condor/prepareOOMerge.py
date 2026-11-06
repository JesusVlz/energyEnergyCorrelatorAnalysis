#!/usr/bin/env python3
"""Audit a complete reserved-donor OO GenGen production and prepare one merge job."""
import argparse
import csv
import datetime as dt
import json
import math
from pathlib import Path
import re
import shutil
from urllib.parse import urlsplit
from zoneinfo import ZoneInfo

from merge_closure_worker import SUFFIXES, require, sha256, validate_manifest, validate_record

CMSSW = '/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_0_11'


def eos_url(path):
    require(str(path).startswith('/eos/'), 'Expected EOS path')
    return 'root://eosuser.cern.ch//' + str(path).lstrip('/')


def local_url(url):
    parsed = urlsplit(url)
    require(parsed.scheme == 'root' and parsed.netloc == 'eosuser.cern.ch', 'Expected CERN EOS user URL')
    path = Path('/' + parsed.path.lstrip('/'))
    require(str(path).startswith('/eos/'), 'Expected EOS URL path')
    return path


def audit_production(production):
    production = Path(production)
    manifest = json.loads((production / 'manifest.json').read_text())
    config = manifest['run_config']
    require(config['sample'] == 'OOMC' and config['modes'] == ['GenGen']
            and config['background_methods'] == 6 and config['data_type'] == 10, 'Unsupported production configuration')
    pool = json.loads((production / 'mixing_pool/manifest.json').read_text())
    provenance = json.loads((production / 'campaign_provenance.json').read_text())
    contract = {'mode': 'GenGen', 'data_type': 10,
                'archive_sha256': manifest['provenance']['archive_sha256'],
                'frozen_card_sha256': sha256(production / 'cards/GenGen.input'),
                'pool_manifest_sha256': sha256(production / 'mixing_pool/manifest.json'),
                'campaign_provenance_sha256': sha256(production / 'campaign_provenance.json'),
                'executable_sha256': provenance['executable_sha256'],
                'signal_headers_sha256': pool['signal_headers_sha256'],
                'source_audit_sha256': pool['source_audit_sha256'],
                'partition_policy': 1, 'seed': 20261030, 'pool_cap': 100000,
                'holdout_policy': 'fixed_reserved_donors_exclude_overlapping_signals',
                'campaign_provenance': provenance}
    require(sha256(production / 'analysis_payload.tar.gz') == contract['archive_sha256'], 'Frozen payload checksum mismatch')
    for file, key in [('signal_headers.manifest.json', 'signal_headers_sha256'),
                      ('source_audit.manifest.json', 'source_audit_sha256')]:
        require(sha256(production / 'mixing_pool' / file) == contract[key], 'Frozen audit checksum mismatch')
    headers = json.loads((production / 'mixing_pool/signal_headers.manifest.json').read_text())
    sources = {row['path']: row for row in headers['files']}
    require(len(sources) == len(headers['files']), 'Duplicate source headers')
    original_inputs = manifest['inputs']
    count = manifest['input_count']
    require(len(original_inputs) == count == manifest['job_count'] > 0, 'Production input/job count mismatch')
    require(len({row['url'] for row in original_inputs}) == count, 'Duplicate production source URL')
    by_url = {row['url']: row for row in original_inputs}
    rows = list(csv.reader((production / 'jobs.tsv').read_text().splitlines(), delimiter='\t'))
    require(len(rows) == count and all(len(row) == 5 for row in rows), 'Incomplete production queue')
    source_dir = Path(config['output_dir']) / 'GenGen'
    entries = []
    totals = dict(all_events=0, inclusive_jet_entries=0, eec_entries=0)
    ledger_totals = {}
    artifact_bytes = 0
    root_bytes = 0
    for mode, input_url, output_url, stem, memory in rows:
        require(mode == 'GenGen' and input_url in by_url, 'Unexpected production queue input')
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', stem) is not None, 'Unsafe production output stem')
        require(output_url == eos_url(source_dir / (stem + '.root')), 'Production output path mismatch')
        source_path = by_url[input_url]['path']
        require(source_path in sources, 'Missing frozen source headers')
        entry = {'stem': stem, 'input_url': input_url, 'output_url': output_url, 'source': sources[source_path]}
        record_path = source_dir / (stem + '.run.json')
        record = json.loads(record_path.read_text())
        validate_record(record, entry, contract)
        entry.update(run_record=record, run_metadata={'sha256': sha256(record_path), 'size_bytes': record_path.stat().st_size})
        for name, artifact in record['artifacts'].items():
            path = source_dir / name
            require(path.stat().st_size == artifact['size_bytes'], 'Published artifact size mismatch: ' + str(path))
            artifact_bytes += artifact['size_bytes']
            if path.suffix == '.root':
                root_bytes += artifact['size_bytes']
        for key in totals:
            value = record['root_validation'][key]
            require(math.isfinite(value) and value >= 0, 'Invalid ROOT exposure')
            totals[key] += value
        for key in ('source_entries', 'selected_signal_events', 'matched_signal_events', 'matched_jets',
                    'contribution_rows', 'reserved_signal_failures', 'partner_assignments'):
            ledger_totals[key] = ledger_totals.get(key, 0) + record['ledger_validation'][key]
        entries.append(entry)
    require({path.name for path in source_dir.glob('*.root')} == {entry['stem'] + '.root' for entry in entries}, 'Unexpected/missing ROOT output')
    require({path.name for path in source_dir.glob('*.run.json')} == {entry['stem'] + '.run.json' for entry in entries}, 'Unexpected/missing run record')
    result = {'schema_version': 1, 'kind': 'OOMC_GenGen_closure_merge',
              'tag': config['tag'], 'production_dir': str(production),
              'source_output_url': eos_url(source_dir), 'input_count': count,
              'contract': contract, 'entries': entries, 'expected_root_totals': totals,
              'ledger_totals': ledger_totals, 'root_bytes': root_bytes, 'artifact_bytes': artifact_bytes,
              'production_manifest_sha256': sha256(production / 'manifest.json'),
              'production_jobs_sha256': sha256(production / 'jobs.tsv'),
              'created_at': dt.datetime.now(ZoneInfo('America/Chicago')).isoformat(),
              'limitations': ['Histogram merging does not establish physics closure or final covariance.',
                              'Split and integrated campaigns repeat physical signals and must remain separate.',
                              'Fixed donors are reused across jobs; original per-job identity and contribution ledgers are retained.',
                              'Preparation verifies success/provenance/size; workers verify every artifact SHA256 before merging.']}
    validate_manifest(result)
    return result


def prepare(production, submit_dir=None, output_dir=None):
    production = Path(production).absolute()
    audited = audit_production(production)
    tag = audited['tag']
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', tag) is not None, 'Unsafe campaign tag')
    workspace = Path(submit_dir) if submit_dir else Path('/afs/cern.ch/user/j/jvelazqu/condor/merge/OOMC') / tag
    destination = Path(output_dir) if output_dir else local_url(audited['source_output_url']).parent / 'merged'
    for path in (workspace, destination):
        require(not any(char.isspace() or char in "\"',;" for char in str(path)), 'Unsupported submission path')
    require(not workspace.exists(), 'Merge campaign already exists: ' + str(workspace))
    require(not destination.exists() or not any(destination.iterdir()), 'Merge output directory is not empty')
    audited.update(destination_url=eos_url(destination), submit_dir=str(workspace))
    source = Path(__file__).parent
    helpers = ('merge_closure_worker.py', 'merge_root.py')
    for name in helpers:
        require((source / name).is_file(), 'Missing merge helper: ' + name)
    audited['helper_sha256'] = {name: sha256(source / name) for name in helpers}
    # Keep enough scratch even without compression; CSVs are removed after entering the archive.
    disk_mb = max(16000, math.ceil((audited['artifact_bytes'] + 2 * audited['root_bytes']) / 1000000 * 1.3) + 2000)
    audited['resources'] = {'cpus': 1, 'memory_mb': 8000, 'disk_mb': disk_mb, 'job_flavour': 'workday'}
    workspace.mkdir(parents=True)
    (workspace / 'logs').mkdir()
    destination.mkdir(parents=True, exist_ok=True)
    for name in helpers:
        shutil.copy2(source / name, workspace / name)
    shutil.copy2(source / 'prepareOOMerge.py', workspace / 'prepareOOMerge.py')
    (workspace / 'manifest.json').write_text(json.dumps(audited, indent=2, allow_nan=False) + '\n')
    manifest_hash = sha256(workspace / 'manifest.json')
    wrapper = f'''#!/bin/bash
set -Eeuo pipefail
trap 'eec_status=$?; echo "MERGE_JOB_FAILED line=$LINENO status=$eec_status" >&2; exit "$eec_status"' ERR
eec_scratch=$(pwd -P)
case "$eec_scratch" in /afs/*|/eos/*) echo 'Worker requires local scratch' >&2; exit 2 ;; esac
echo '{manifest_hash}  manifest.json' | sha256sum -c -
'''
    for name, digest in audited['helper_sha256'].items():
        wrapper += f"echo '{digest}  {name}' | sha256sum -c -\n"
    wrapper += f'''set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cd {CMSSW}
eec_runtime=$(scram runtime -sh)
eval "$eec_runtime"
set -u
cd "$eec_scratch"
python3 merge_closure_worker.py --manifest manifest.json --work-dir "$eec_scratch/merge-work"
'''
    (workspace / 'run_merge.sh').write_text(wrapper)
    (workspace / 'run_merge.sh').chmod(0o755)
    transfers = ','.join(str(workspace / name) for name in ('manifest.json', *helpers))
    (workspace / 'jobs.sub').write_text(f'''universe = vanilla
initialdir = {workspace}
executable = {workspace}/run_merge.sh
log = {workspace}/logs/jobs.$(ClusterId).log
output = {workspace}/logs/GenGen.$(ClusterId).$(ProcId).out
error = {workspace}/logs/GenGen.$(ClusterId).$(ProcId).err
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {transfers}
transfer_output_files = ""
request_cpus = 1
request_memory = 8000MB
request_disk = {disk_mb}MB
requirements = (OpSysAndVer =?= "AlmaLinux9") && (Arch =?= "X86_64")
+JobFlavour = "workday"
+JobBatchName = "OOMC_merge_{tag}"
queue 1
''')
    # No blind retries: publication refuses to replace an existing artifact.
    (workspace / 'submit.sh').write_text(f'#!/bin/bash\nset -euo pipefail\ncd {workspace}\nexec condor_submit "$@" jobs.sub\n')
    (workspace / 'submit.sh').chmod(0o755)
    (workspace / 'README.txt').write_text(
        f'OO GenGen merge: {tag}\nProduction: {production}\nInputs: {audited["input_count"]}\n'
        f'Output directory: {destination}\n\n'
        'The worker checks every run record and five artifact SHA256 values, then merges every ROOT input.\n'
        'There is no skip-bad-files option. Histogram sums, axes and invariant JCard settings are checked.\n'
        'Original per-job CSVs and run records are retained in merged_OOMC_GenGen.ledgers.tar.gz.\n'
        'Consume merged_OOMC_GenGen.root only with the final successful merged_OOMC_GenGen.run.json.\n'
        'The report records SHA256 and size for the ROOT and ledger archive. Existing outputs are protected.\n'
        'Physics closure and donor/cross-job covariance remain open. Keep split/integrated separate.\n\n'
        f'Dry run: cd {workspace} && ./submit.sh -dry-run dryrun.ads\n'
        f'Submit: cd {workspace} && ./submit.sh\n')
    print(json.dumps({'campaign': str(workspace), 'output_dir': str(destination),
                      'inputs': audited['input_count'], 'expected_root_totals': audited['expected_root_totals'],
                      'resources': audited['resources']}), flush=True)
    return audited


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--production', required=True, type=Path)
    parser.add_argument('--submit-dir', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--audit-only', action='store_true')
    args = parser.parse_args()
    if args.audit_only:
        result = audit_production(args.production)
        print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, indent=2))
    else:
        prepare(args.production, args.submit_dir, args.output_dir)


if __name__ == '__main__':
    main()
