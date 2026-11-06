#!/usr/bin/env python3
"""Validate complete independent folds, sum raw response counts, then project.

No normalization, unfolding or iteration selection is performed. The manifest
gate and physical-event uniqueness audit precede every histogram merge.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile

from response_worker import REPO, card_values, digest, event_key, latest_objects, raw_report


def canonical(path):
    return str(path).replace('/eos/user/j/jvelazqu/', '/eos/home-j/jvelazqu/')


def validate_inputs(folders, system, manifest, database, allow_bounded=False):
    tasks = {canonical(t['input']): t for t in manifest['tasks'] if t['system'] == system} if manifest else None
    if tasks is None and not allow_bounded:
        raise ValueError('Production merge requires the full frozen input manifest')
    connection = sqlite3.connect(str(database))
    connection.execute('CREATE TABLE identities (namespace TEXT, run TEXT, lumi TEXT, event TEXT, PRIMARY KEY(namespace,run,lumi,event)) WITHOUT ROWID')
    records, seen_inputs, expected_card, common_hashes = [], set(), None, None
    try:
        for folder in folders:
            receipt_path = folder/'receipt.json'
            record = json.loads(receipt_path.read_text())
            if record['status'] != 'COMPLETE_ENGINEERING_ONLY' or record['system'] != system:
                raise ValueError('Incomplete or wrong-system response worker: '+str(folder))
            skim = record['skim']
            source = canonical(record['source_input'])
            if source in seen_inputs:
                raise ValueError('Duplicate response source input: '+source)
            seen_inputs.add(source)
            if not allow_bounded and skim['processed_entries'] != skim['input_entries']:
                raise ValueError('Truncated response source input: '+source)
            if tasks is not None:
                if source not in tasks:
                    raise ValueError('Response input is absent from frozen manifest: '+source)
                task = tasks[source]
                identity = task['input_identity']
                if (skim['input_entries'] != identity['entries'] or skim['processed_entries'] != identity['entries']
                        or skim['input_size_bytes'] != identity['size_bytes'] or skim['sample_namespace'] != task['sample_namespace']):
                    raise ValueError('Response input coverage/identity differs from frozen manifest: '+source)
            hashes = {name: record[name] for name in ('native_sha256', 'card_sha256', 'helper_sha256', 'pair_counts_sha256')}
            if common_hashes is not None and hashes != common_hashes:
                raise ValueError('Mixed response source/binary/card versions')
            common_hashes = hashes
            ledger = folder/'event_identity.csv.gz'
            if digest(ledger) != skim['ledger_sha256']:
                raise ValueError('Changed event-identity ledger: '+str(folder))
            counts = [0, 0]
            with gzip.open(ledger, 'rt') as stream:
                for entry, row in enumerate(csv.DictReader(stream)):
                    identity = tuple(int(row[k]) for k in ('run', 'lumi', 'event'))
                    key = event_key(dict(sample=system, sample_namespace=skim['sample_namespace'], id=identity))
                    which = key % 2
                    if (int(row['source_entry']) != entry or int(row['event_key']) != key or
                            int(row['fold']) != which or int(row['fold_entry']) != counts[which]):
                        raise ValueError('Incorrect physical-event fold ledger: '+str(folder))
                    try:
                        connection.execute('INSERT INTO identities VALUES(?,?,?,?)', (skim['sample_namespace'], *map(str, identity)))
                    except sqlite3.IntegrityError as error:
                        raise ValueError('Duplicate global physical event: '+repr(identity)) from error
                    counts[which] += 1
            connection.commit()
            if counts != skim['fold_entries'] or sum(counts) != skim['processed_entries']:
                raise ValueError('Incomplete fold ledger: '+str(folder))
            for label in ('train', 'validation'):
                path = folder/(label+'_raw.root')
                info = record['folds'][label]
                if digest(path) != info['sha256']:
                    raise ValueError('Changed native response output: '+str(path))
                actual = raw_report(path)
                if expected_card is not None and actual['card'] != expected_card:
                    raise ValueError('Incompatible raw response JCards')
                expected_card = actual['card']
                if not actual['inclusive_activity_accounted']:
                    raise ValueError('Inclusive activity flow is present; retain/account it before projection, without an activity cut')
            records.append(dict(folder=str(folder), receipt_sha256=digest(receipt_path), source_input=source,
                                entries=skim['processed_entries'], fold_entries=counts, hashes=hashes))
        if tasks is not None and seen_inputs != set(tasks):
            raise ValueError('Incomplete frozen source coverage; missing '+str(len(set(tasks)-seen_inputs))+' files')
        if not records:
            raise ValueError('No response inputs')
        return dict(system=system, inputs=records, files=len(records),
                    events=sum(r['entries'] for r in records),
                    fold_events=[sum(r['fold_entries'][i] for r in records) for i in (0, 1)],
                    source_sha256=common_hashes, independent_identity_audit='PASS: exact namespace/run/lumi/event unique across all input ledgers',
                    bounded_verification_only=allow_bounded)
    finally:
        connection.close()


def histogram_objects(directory, prefix=''):
    """Visit the latest histogram cycles, retaining their directory paths."""
    for name, obj in latest_objects(directory):
        if not prefix and name == 'JCard':
            continue  # Metadata is copied once, never summed.
        path = prefix + name
        if obj.InheritsFrom('TDirectory'):
            yield from histogram_objects(obj, path + '/')
        elif obj.InheritsFrom('TH1') or obj.InheritsFrom('THnBase'):
            yield path, obj


def histogram_layout(histogram):
    """Require identical axes instead of letting ROOT rebin or relabel inputs."""
    if histogram.InheritsFrom('THnBase'):
        axes = [histogram.GetAxis(i) for i in range(histogram.GetNdimensions())]
    else:
        axes = [histogram.GetXaxis(), histogram.GetYaxis(), histogram.GetZaxis()][:histogram.GetDimension()]
    return (histogram.ClassName(), tuple(
        (tuple(axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)),
         tuple(axis.GetBinLabel(i) for i in range(1, axis.GetNbins() + 1)))
        for axis in axes))


def merge_raw(paths, output, provenance):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if not paths:
        raise ValueError('No raw response inputs')
    accumulators, layouts = {}, {}
    first = ROOT.TFile.Open(str(paths[0]), 'READ')
    if not first or first.IsZombie():
        raise ValueError('Unreadable first raw response input')
    try:
        original_card = card_values(first.GetDirectory('JCard'))
        for name, obj in histogram_objects(first):
            clone = obj.Clone(obj.GetName())
            if clone.InheritsFrom('TH1'):
                clone.SetDirectory(0)
            accumulators[name] = clone
            layouts[name] = histogram_layout(obj)
        if not accumulators:
            raise ValueError('No histograms in raw response input')
        for path in paths[1:]:
            handle = ROOT.TFile.Open(str(path), 'READ')
            if not handle or handle.IsZombie():
                raise ValueError('Unreadable raw response input '+str(path))
            try:
                if card_values(handle.GetDirectory('JCard')) != original_card:
                    raise ValueError('Incompatible raw response card during merge')
                incoming = dict(histogram_objects(handle))
                if set(incoming) != set(accumulators):
                    raise ValueError('Different raw histogram paths in '+str(path)+
                                     ': missing='+repr(sorted(set(accumulators)-set(incoming)))+
                                     ', unexpected='+repr(sorted(set(incoming)-set(accumulators))))
                for name, histogram in accumulators.items():
                    if histogram_layout(incoming[name]) != layouts[name]:
                        raise ValueError('Incompatible raw histogram class/binning/labels: '+name)
                    added = histogram.Add(incoming[name])
                    # TH1::Add returns bool; THnBase::Add returns void. The
                    # latter's class and axis checks are enforced above.
                    if added is not None and not added:
                        raise ValueError('ROOT failed to add raw histogram '+name)
            finally:
                handle.Close()
        target = ROOT.TFile.Open(str(output), 'RECREATE')
        if not target or target.IsZombie():
            raise ValueError('Cannot create merged raw response output')
        try:
            for path, histogram in accumulators.items():
                parts = path.split('/')
                directory = target
                for name in parts[:-1]:
                    child = directory.GetDirectory(name)
                    directory = child if child else directory.mkdir(name)
                directory.cd()
                if histogram.Write(parts[-1]) <= 0:
                    raise ValueError('Cannot write merged raw histogram '+path)
            card = target.mkdir('JCard')
            for name, obj in latest_objects(first.GetDirectory('JCard')):
                card.cd()
                obj.Write(name)
            target.cd()
            ROOT.TObjString(json.dumps(provenance, sort_keys=True)).Write('ResponseMergeProvenance')
        finally:
            target.Close()
    finally:
        first.Close()
    return raw_report(output)


def inclusive_input(original, output):
    """Request an exact union of existing activity bins; histograms stay intact."""
    import ROOT
    shutil.copyfile(original, output)
    handle = ROOT.TFile.Open(str(output), 'UPDATE')
    try:
        card = handle.GetDirectory('JCard')
        values = card_values(card)['CentralityBinEdges']
        if values != [-1., 0., 30., 50., 90., 100.]:
            raise ValueError('Unexpected OO activity edges; cannot form exact inclusive union')
        native_card = handle.mkdir('SourceJCard')
        for name, obj in latest_objects(card):
            native_card.cd()
            obj.Write(name)
        card.cd()
        vector = ROOT.TVectorF(1, 2)
        vector.GetMatrixArray()[0] = values[0]
        vector.GetMatrixArray()[1] = values[-1]
        vector.Write('CentralityBinEdges', ROOT.TObject.kOverwrite)
        handle.cd()
        ROOT.TObjString(json.dumps(dict(activity='Inclusive all otherwise-selected events; no activity event cut',
                                       original_edges=values, projected_edges=[values[0], values[-1]],
                                       prescription='Exact union of existing histogram bins; activity flow verified empty before union'))).Write('ActivityUnionProvenance')
    finally:
        handle.Close()


def project_with_metadata(raw, output, power, scratch, source_raw):
    command = [sys.executable, str(REPO/'validation/run3Unfolding/project_response.py'), str(raw), str(output),
               '--power', str(power), '--mask', '108803', '--scratch', str(scratch)]
    with output.with_suffix('.log').open('w') as log:
        subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, check=True)
    import ROOT
    target = ROOT.TFile.Open(str(output), 'UPDATE')
    source = ROOT.TFile.Open(str(source_raw), 'READ')
    try:
        # InputJCard already retains the exact card passed to the unchanged
        # projection macro. SourceJCard preserves the native activity partition.
        directory = target.mkdir('SourceJCard')
        for name, obj in latest_objects(source.GetDirectory('JCard')):
            directory.cd()
            obj.Write(name)
        target.cd()
        source.Get('ResponseMergeProvenance').Write('ResponseMergeProvenance')
    finally:
        target.Close()
        source.Close()
    check = ROOT.TFile.Open(str(output), 'READ')
    try:
        if not check or check.IsZombie() or check.TestBit(ROOT.TFile.kRecovered) or not check.GetDirectory('SourceJCard'):
            raise ValueError('Invalid final projection metadata')
    finally:
        if check:
            check.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path, help='JSON list of local worker output directories')
    parser.add_argument('output', type=Path, help='New merged/projected result directory')
    parser.add_argument('--system', choices=('OO', 'pp'), required=True)
    parser.add_argument('--manifest', type=Path, help='Full frozen all-MC input manifest; required in production')
    parser.add_argument('--allow-bounded', action='store_true', help='Explicit engineering fixture/smoke mode only')
    parser.add_argument('--no-project', action='store_true', help='Validated raw matrices only')
    parser.add_argument('--scratch-parent', type=Path)
    args = parser.parse_args()
    folders = [Path(p).resolve(strict=True) for p in json.loads(args.inputs.read_text())]
    manifest = json.loads(args.manifest.read_text()) if args.manifest else None
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='eec-response-merge-', dir=args.scratch_parent) as directory:
        scratch = Path(directory)
        report = validate_inputs(folders, args.system, manifest, scratch/'identities.sqlite', args.allow_bounded)
        report['manifest_sha256'] = digest(args.manifest) if args.manifest else None
        report['normalization'] = 'Raw weighted counts summed before projection; inherited generator unity policy is unvalidated'
        report['outputs'] = {}
        for label in ('train', 'validation'):
            raw = scratch/(args.system+'_'+label+'_raw.root')
            raw_qa = merge_raw([p/(label+'_raw.root') for p in folders], raw, dict(report, fold=label))
            outputs = [raw]
            if not args.no_project:
                inputs = [('split' if args.system=='OO' else 'inclusive', raw)]
                if args.system == 'OO':
                    inclusive = scratch/(args.system+'_'+label+'_inclusive_input.root')
                    inclusive_input(raw, inclusive)
                    inputs.append(('inclusive', inclusive))
                for activity, input_path in inputs:
                    for power in (1, 2):
                        projected = scratch/(args.system+'_'+label+'_'+activity+'_N'+str(power)+'.root')
                        project_with_metadata(input_path, projected, power, scratch/(projected.stem+'_build'), raw)
                        outputs.append(projected)
            for path in outputs:
                target = output/path.name
                checksum = digest(path)
                shutil.copyfile(path, target)
                if digest(target) != checksum:
                    raise ValueError('Merged output publication checksum mismatch')
                report['outputs'][path.name] = dict(sha256=checksum, bytes=target.stat().st_size)
            report[label+'_qa'] = raw_qa
        report['status'] = 'COMPLETE_RESPONSE_MATRICES_ENGINEERING_ONLY'
        (output/'merge_receipt.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], system=args.system, files=report['files'], events=report['events'])))


if __name__ == '__main__':
    main()
