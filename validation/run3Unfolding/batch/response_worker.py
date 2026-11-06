#!/usr/bin/env python3
"""Full-file native response matrices on the existing independent event folds.

Only synchronized forest skimming and orchestration are added. Native selections,
weights and response construction remain unchanged. This never unfolds.
"""
import argparse
from array import array
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO/'validation/run3Pairs'))
from pair_counts import event_key, fold

NAMESPACE = {
    'OO': 'OO2025_PHOO_pythia_embedded_2GeV_constit_cut_260926',
    'pp': 'ppRef2024_Pythia_QCD_2GeV_constit_cut_260926',
}
RESPONSE_NAMES = ('jetPtUnfoldingMeasured', 'jetPtUnfoldingTruth',
    'jetPtUnfoldingResponse', 'oneDimensionalJetPtUnfoldingMeasured',
    'oneDimensionalJetPtUnfoldingTruth', 'oneDimensionalJetPtUnfoldingResponse',
    'jetPtUnfoldingCovariance')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def card_values(directory):
    if not directory:
        raise ValueError('Missing embedded JCard')
    result = {}
    for key in directory.GetListOfKeys():
        name = key.GetName()
        if name in result:
            continue
        obj = directory.Get(name)
        if obj.ClassName().startswith('TVectorT'):
            data = obj.GetMatrixArray()
            result[name] = [float(data[i]) for i in range(obj.GetNrows())]
        elif obj.InheritsFrom('TObjString'):
            result[name] = str(obj.GetString())
    return result


def latest_objects(directory):
    # HiForest contains older partial TTree cycles. Native readers use Get(name),
    # so explicitly select the same latest cycle, never append earlier cycles.
    seen = set()
    for key in directory.GetListOfKeys():
        name = key.GetName()
        if name not in seen:
            seen.add(name)
            yield name, directory.Get(name)


def skim(input_path, output_dir, system, namespace, max_events=0, source_input=None):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    source = ROOT.TFile.Open(str(input_path), 'READ')
    if not source or source.IsZombie() or source.TestBit(ROOT.TFile.kRecovered):
        raise ValueError('Unreadable/recovered input forest: '+str(input_path))
    event_tree = source.Get('hiEvtAnalyzer/HiTree')
    if not event_tree:
        raise ValueError('Missing hiEvtAnalyzer/HiTree')
    total = int(event_tree.GetEntries())
    count = min(total, max_events) if max_events else total
    if count <= 0 or not namespace:
        raise ValueError('Empty input or missing physical sample namespace')
    for name in ('run', 'lumi', 'evt'):
        if not event_tree.GetBranch(name):
            raise ValueError('Missing physical identity branch: '+name)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [output_dir/(label+'_forest.root') for label in ('train', 'validation')]
    if any(p.exists() for p in paths):
        raise ValueError('Skim outputs must be new')
    assignments = ROOT.std.vector('unsigned char')()
    assignments.reserve(count)
    selected = [0, 0]
    identities = set()
    event_tree.SetBranchStatus('*', 0)
    for name in ('run', 'lumi', 'evt'):
        event_tree.SetBranchStatus(name, 1)
    ledger = output_dir/'event_identity.csv.gz'
    with gzip.open(ledger, 'wt') as stream:
        stream.write('source_entry,run,lumi,event,event_key,fold,fold_entry\n')
        for entry in range(count):
            if event_tree.GetEntry(entry) <= 0:
                raise ValueError('Unreadable event identity at entry '+str(entry))
            identity = (int(event_tree.run), int(event_tree.lumi), int(event_tree.evt))
            if identity in identities or any(v < 0 or v >= 1 << 64 for v in identity):
                raise ValueError('Duplicate/invalid physical identity: '+repr(identity))
            identities.add(identity)
            event = dict(sample=system, sample_namespace=namespace, id=identity)
            which = fold(event)
            assignments.push_back(which)
            stream.write(','.join(map(str, (entry, *identity, event_key(event), which, selected[which])))+'\n')
            selected[which] += 1
    event_tree.SetBranchStatus('*', 1)
    if not hasattr(ROOT, 'SplitResponseTree'):
        ROOT.gInterpreter.Declare('''
        #include <vector>
        #include <stdexcept>
        #include "TTree.h"
        void SplitResponseTree(TTree* input, TTree* train, TTree* validation,
                               const std::vector<unsigned char>& folds) {
          for (Long64_t entry=0; entry<(Long64_t)folds.size(); ++entry) {
            if (input->GetEntry(entry)<0) throw std::runtime_error("Unreadable forest entry");
            TTree* output=folds[entry] ? validation : train;
            if(output->Fill()<0) throw std::runtime_error("Cannot fill synchronized skim");
          }
        }
        ''')
    outputs = [ROOT.TFile.Open(str(p), 'RECREATE') for p in paths]
    records = []
    try:
        def copy_directory(indir, outdirs, prefix=''):
            for name, obj in latest_objects(indir):
                path = prefix+name
                if obj.InheritsFrom('TDirectory'):
                    copy_directory(obj, [out.mkdir(name) for out in outdirs], path+'/')
                elif obj.InheritsFrom('TTree'):
                    entries = int(obj.GetEntries())
                    metadata = path.startswith('HiForestInfo/')
                    if not metadata and entries != total:
                        raise ValueError('Unsynchronized event tree '+path+': '+str(entries)+' != '+str(total))
                    clones = []
                    for out in outdirs:
                        out.cd()
                        clone = obj.CloneTree(-1 if metadata else 0, 'fast' if metadata else '')
                        clone.SetDirectory(out)
                        clones.append(clone)
                    if not metadata:
                        ROOT.SplitResponseTree(obj, *clones, assignments)
                    for out, clone in zip(outdirs, clones):
                        out.cd()
                        clone.Write(name, ROOT.TObject.kOverwrite)
                    records.append(dict(path=path, input_entries=entries,
                                        output_entries=[int(c.GetEntries()) for c in clones], metadata=metadata))
                else:
                    for out in outdirs:
                        out.cd()
                        obj.Write(name)
        copy_directory(source, outputs)
    finally:
        for output in outputs:
            output.Close()
        source.Close()
    # Reopen closed files; metadata and each synchronized tree must survive.
    for which, path in enumerate(paths):
        handle = ROOT.TFile.Open(str(path), 'READ')
        if not handle or handle.IsZombie() or handle.TestBit(ROOT.TFile.kRecovered):
            raise ValueError('Unreadable/recovered synchronized skim')
        try:
            for record in records:
                tree = handle.Get(record['path'])
                if not tree or int(tree.GetEntries()) != record['output_entries'][which]:
                    raise ValueError('Skim tree count changed on reopening: '+record['path'])
        finally:
            handle.Close()
    return dict(input=str(source_input or input_path), staged_input=str(input_path),
                input_size_bytes=Path(input_path).stat().st_size,
                input_entries=total, processed_entries=count, fold_entries=selected,
                sample_namespace=namespace, fold_algorithm='pair_counts.fold: namespaced FNV1a64 modulo2; 0 train,1 validation',
                ledger_sha256=digest(ledger), trees=records,
                skim_bytes=[p.stat().st_size for p in paths]), paths


def raw_report(path):
    import ROOT
    handle = ROOT.TFile.Open(str(path), 'READ')
    if not handle or handle.IsZombie() or handle.TestBit(ROOT.TFile.kRecovered):
        raise ValueError('Unreadable/recovered native response output')
    histograms = {}
    try:
        actual = card_values(handle.GetDirectory('JCard'))
        required = dict(WeightExponent=[1., 2.], MinJetPtUnfoldingTruth=[5.],
                        FirstUnfoldedJetPtBin=[3.], LastUnfoldedJetPtBin=[5.],
                        MaxEventsPerFile=[0.], MinTrackPtCut=[2.], McCorrelationType=[1.])
        for key, expected in required.items():
            if actual.get(key) != expected:
                raise ValueError('Unexpected embedded response selection '+key)
        for name in RESPONSE_NAMES:
            histogram = handle.Get(name)
            if not histogram or not histogram.InheritsFrom('THnSparse'):
                raise ValueError('Missing response output '+name)
            dimensions = histogram.GetNdimensions()
            coordinates = array('i', [0]*dimensions)
            activity_axis = ((2 if name.endswith('Response') else 1) if name.startswith('oneDimensional')
                             else (3 if name.endswith(('Response', 'Covariance')) else 2))
            report = dict(entries=histogram.GetEntries(), sumw_all=0., sumw_activity_flow=0.,
                          activity_flow_cells=0, negative_cells=0, nonfinite_cells=0)
            for index in range(histogram.GetNbins()):
                value = float(histogram.GetBinContent(index, coordinates))
                if not math.isfinite(value):
                    report['nonfinite_cells'] += 1
                    continue
                report['negative_cells'] += int(value < 0)
                report['sumw_all'] += value
                if (coordinates[activity_axis] < 1 or coordinates[activity_axis] > histogram.GetAxis(activity_axis).GetNbins()) and value:
                    report['activity_flow_cells'] += 1
                    report['sumw_activity_flow'] += value
            if report['nonfinite_cells']:
                raise ValueError('Nonfinite response weight in '+name)
            histograms[name] = report
        return dict(card=actual, histograms=histograms,
                    inclusive_activity_accounted=all(not h['activity_flow_cells'] for h in histograms.values()),
                    normalization='Inherited generator factor1; cross-section normalization unvalidated',
                    covariance='Established within-jet second moments only; full event covariance unvalidated',
                    detector_calibration='Pure generator-particle response; no detector correction applied')
    finally:
        handle.Close()


def run(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    card = Path(__file__).parent/'cards'/(args.system+'_R04_response.input')
    analysis_repo = (args.analysis_repo or REPO).resolve(strict=True)
    native = analysis_repo/'eecAnalysis'
    started = time.time()
    receipt = dict(status='RUNNING', system=args.system, source_input=str(args.source_input or args.input),
                   card=str(card), card_sha256=digest(card), native_sha256=digest(native),
                   helper_sha256=digest(__file__), pair_counts_sha256=digest(REPO/'validation/run3Pairs/pair_counts.py'))
    record_path = output/'receipt.json'
    record_path.write_text(json.dumps(receipt, indent=2)+'\n')
    try:
        with tempfile.TemporaryDirectory(prefix='eec-response-worker-', dir=args.scratch_parent) as scratch:
            work = Path(scratch)
            provenance, forests = skim(args.input.resolve(strict=True), work, args.system,
                                      args.namespace or NAMESPACE[args.system], args.max_events, args.source_input)
            receipt['skim'] = provenance
            shutil.copyfile(work/'event_identity.csv.gz', output/'event_identity.csv.gz')
            receipt['folds'] = {}
            for index, label in enumerate(('train', 'validation')):
                if provenance['fold_entries'][index] <= 0:
                    raise ValueError('Empty physical-event fold')
                file_list = work/(label+'.txt')
                file_list.write_text(str(forests[index])+'\n')
                raw = work/(label+'_raw.root')
                command = [str(native), str(file_list), str(card), str(raw), '0', 'true']
                with (output/(label+'.log')).open('w') as log:
                    subprocess.run(command, cwd=analysis_repo, stdout=log, stderr=subprocess.STDOUT, check=True)
                report = raw_report(raw)
                report['command'] = command
                metadata = dict(receipt, skim=provenance, fold=index, label=label,
                                scope='Response matrices only; no unfolding or iteration choice')
                import ROOT
                handle = ROOT.TFile.Open(str(raw), 'UPDATE')
                ROOT.TObjString(json.dumps(metadata, sort_keys=True)).Write('ResponseProvenance')
                handle.Close()
                # Verify after the final UPDATE while the file is still local.
                raw_report(raw)
                target = output/raw.name
                expected = digest(raw)
                shutil.copyfile(raw, target)
                if digest(target) != expected:
                    raise ValueError('Response output publication checksum mismatch')
                report.update(output=str(target), sha256=expected, bytes=target.stat().st_size)
                receipt['folds'][label] = report
            receipt.update(status='COMPLETE_ENGINEERING_ONLY', elapsed_seconds=time.time()-started)
    except Exception as error:
        receipt.update(status='FAILED', error=str(error), elapsed_seconds=time.time()-started)
        record_path.write_text(json.dumps(receipt, indent=2)+'\n')
        raise
    record_path.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(dict(status=receipt['status'], output=str(output), fold_entries=receipt['skim']['fold_entries'])), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path, help='New output folder; staging stays local')
    parser.add_argument('--system', choices=('OO', 'pp'), required=True)
    parser.add_argument('--analysis-repo', type=Path, help='Runtime analysis checkout containing executable and calibration files')
    parser.add_argument('--namespace', help='Exact physical sample namespace; defaults to established campaign')
    parser.add_argument('--source-input', type=Path, help='Original input path if INPUT was staged locally')
    parser.add_argument('--max-events', type=int, default=0, help='0=all entries; positive only for bounded verification')
    parser.add_argument('--scratch-parent', type=Path, help='Existing local scratch parent')
    args = parser.parse_args()
    if args.max_events < 0:
        parser.error('--max-events must be nonnegative')
    run(args)


if __name__ == '__main__':
    main()
