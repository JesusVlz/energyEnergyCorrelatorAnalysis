#!/usr/bin/env python3
"""Read ROOT headers and metadata for existing OO GenGen campaign inputs; no event scan."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys
import time

from pilot_input_audit import inspect_file, digest, PRODUCER_ROOT, PRODUCER_FILES

GEN_JETS = ('ngen', 'genpt', 'geneta', 'genphi', 'WTAgeneta', 'WTAgenphi')
PATTERN = re.compile(r'unembedded_L1HLT_objects_MC_([0-9]+)\.root')


def audit(directory, output, expected_files=333):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    files = sorted((p for p in directory.iterdir() if PATTERN.fullmatch(p.name) and p.is_file()),
                   key=lambda p: int(PATTERN.fullmatch(p.name)[1]))
    records, schema_groups, failures = [], {}, []
    if len(files) != expected_files:
        failures.append({'problem': 'Discovered file count differs from requested existing campaign',
                         'expected': expected_files, 'actual': len(files)})
    start = time.monotonic()
    for index, path in enumerate(files):
        try:
            metadata = inspect_file(ROOT, str(path))
            schema = metadata.pop('schema')
            missing = sorted(set(GEN_JETS)-set(schema['akCs4PFJetAnalyzer/t']))
            if 'pthat' not in schema['hiEvtAnalyzer/HiTree']:
                missing.append('hiEvtAnalyzer/HiTree:pthat')
            if missing:
                raise ValueError('Missing native generator/MC branches: '+str(missing))
            if metadata['source_entries'] <= 0:
                raise ValueError('Empty event trees')
            encoded = json.dumps(schema, sort_keys=True).encode()
            schema_key = hashlib.sha256(encoded).hexdigest()
            if schema_key not in schema_groups:
                schema_groups[schema_key] = {'schema': schema, 'files': []}
            schema_groups[schema_key]['files'].append(str(path))
            metadata.update(status='PASS_REQUIRED_HEADERS', schema_sha256=schema_key,
                            physical_identity='run UInt_t / lumi UInt_t / evt ULong64_t',
                            native_generator_collection='akCs4PFJetAnalyzer/t: genpt + selected WTAgeneta/WTAgenphi',
                            generator_event_rows_read=0,
                            jet_constituent_retention='Unknown from headers; production name may describe a jet-constituent cut and is not proof of stored-particle truncation')
            records.append(metadata)
        except Exception as error:
            failure = {'path': str(path), 'status': 'FAIL', 'problem': str(error)}
            failures.append(failure)
            records.append(failure)
            print('SCHEMA FAILURE '+str(path)+': '+str(error), flush=True)
        if (index+1) % 25 == 0 or index+1 == len(files):
            print('Headers audited '+str(index+1)+'/'+str(len(files)), flush=True)
    producer_files = {}
    for relative in PRODUCER_FILES:
        source = PRODUCER_ROOT/relative
        if source.is_file():
            text = source.read_text()
            producer_files[str(source)] = {'sha256': digest(source),
                'relevant_lines': [line.strip() for line in text.splitlines()
                                   if any(word in line for word in ('stableOnly', 'ptMin', 'etaMax', 'GenParticle', 'sube', '2GeV'))][:60]}
    source_counts = Counter(record.get('source_entries', 0) for record in records)
    result = {'status': 'FAIL_REQUIRED_INPUT_SCHEMA' if failures else 'PASS_REQUIRED_HEADERS_ONLY',
              'input_directory': str(directory.resolve()), 'expected_files': expected_files,
              'discovered_files': len(files), 'passed_files': len(files)-sum(record.get('status')=='FAIL' for record in records),
              'total_source_entries': sum(record.get('source_entries', 0) for record in records),
              'total_source_bytes': sum(record.get('size_bytes', 0) for record in records),
              'source_entries_distribution': dict(source_counts), 'schema_groups': schema_groups,
              'files': records, 'failures': failures, 'producer_files': producer_files,
              'audit_scope': 'Read-only filesystem stats, ROOT tree headers/GetEntries/schema, HiForestInfo metadata first entry. Zero generator/event/jet physics rows scanned.',
              'requirements_checked': ['nonempty/openable ROOT', 'aligned event/gen/skim/native R4 trees', 'full physical run/lumi/64-bit event identity schema',
                  'MC pthat/weight and raw hiHF_pf/vz/hiBin', 'native genpt and WTAgeneta/WTAgenphi', 'native reconstructed PF Jet ID/veto branches',
                  'OO primary-vertex/HF-coincidence/pileup filters', 'generator pt/eta/phi/charge/subevent vectors'],
              'limitations': ['No campaign physical-ID duplicate or background-provenance event scan performed.',
                  'No campaign-wide subevent-label values, physical centrality acceptance, particle pT/eta retention, or reconstructed veto losses inferred from headers.',
                  'Local producer stableOnly=True, ptMin=0.25, etaMax=2.5 and absent sta corroborate the source convention; exact production configuration linkage remains unverified.',
                  'A production dataset name containing 2GeV_constit_cut does not establish loss of low-pT stored generator particles; measured low-pT pilot rows remain the relevant bounded evidence.',
                  'Schema success does not establish physics closure, campaign event independence, final covariance, or permission to freeze the estimator.'],
              'root_version': ROOT.gROOT.GetVersion(), 'python_version': sys.version,
              'tool_sha256': digest(__file__), 'elapsed_seconds': time.monotonic()-start}
    with output.open('x') as destination:
        json.dump(result, destination, indent=2, allow_nan=False)
        destination.write('\n')
    print(json.dumps({key: result[key] for key in ('status','discovered_files','passed_files','total_source_entries','elapsed_seconds')}), flush=True)
    return not failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-files', type=int, default=333)
    args = parser.parse_args()
    raise SystemExit(0 if audit(args.input_dir, args.output, args.expected_files) else 1)


if __name__ == '__main__':
    main()
