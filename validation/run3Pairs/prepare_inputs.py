#!/usr/bin/env python3
"""Prepare, but do not submit, the existing AN7.5.1 Run-3 pair workflow.

No MC locations are assumed. New campaigns use user-provided local/EOS-mounted
file lists, existing native export and mergeable batch reductions. Training
counts are merged globally before normalization; holdout events stay separate.
"""
import argparse
import json
import math
import shlex
import sys
from pathlib import Path
from preparation_contract import REPO, VERSION, VERTEX, payload_path, digest


def settings(system, selection, correlation=0):
    return {'DataType': 10 if system == 'OO' else 12,
            'Run3TrackSelection': selection, 'MaxTrackPtRelativeError': 0.1,
            'VertexMaxDistance': VERTEX[selection], 'TrackEtaCut': 2.4,
            'MinTrackPtCut': 0.5, 'MaxTrackPtCut': 500, 'MinPairTrackPt': 2,
            'McCorrelationType': correlation, 'MatchJets': 0, 'SubeventCut': 2,
            'JetRadius': 0.4, 'JetAxis': 1, 'JetEtaCut': 1.6,
            'MinJetPtCut': 40, 'MaxJetPtCut': 160, 'JetPtBinEdgesEEC': '40 80 120 160',
            'WeightExponent': '1 2', 'BackgroundMethods': 0, 'MegaSkimMixing': 0,
            'HFEnergyShift': 0, 'JetPtWeight': 0, 'DisableTrackPairEfficiencyCorrection': 1,
            'SmearDeltaR': 0, 'SmearEnergyWeight': 0, 'TrackEfficiencyVariation': 0,
            'SkipCovarianceMatrix': 1, 'FilledHistograms': 23, 'MaxEventsPerFile': 10000,
            'CentralityBinEdges': '0 30 50 90 100' if system == 'OO' else '-1 100'}


def write_cards(output):
    from export_forest import card_text
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    result = {}
    for system in ('OO', 'pp'):
        base = REPO / ('cardEECOOMCInclusiveR04.input' if system == 'OO' else 'cardEECpp5v36TeVMCR04.input')
        for selection in VERTEX:
            for population, correlation in (('RecoReco', 0), ('RecoGen', 1)):
                path = output / f'{system}_{selection}_{population}.input'
                text = card_text(base, settings(system, selection, correlation))
                # Never inherit a candidate correction/export path from a base card.
                text = '\n'.join(line for line in text.splitlines() if not line.split('#',1)[0].split() or
                                  line.split()[0] not in ('Run3TrackPairCorrectionFile', 'PairExportFile', 'PairExportNamespace')) + '\n'
                path.write_text('# Pair derivation only: ' + VERSION + '; no production payload enabled.\n' + text)
                result[(system, selection, population)] = path.resolve()
    return result


def axes(hist):
    result = []
    for axis in (hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis())[:hist.GetDimension()]:
        result.append([float(axis.GetBinLowEdge(i)) for i in range(1, axis.GetNbins()+2)])
    return result


def audit_tables():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    report = {'status': 'TECHNICAL_STRUCTURE_CHECK_ONLY', 'physics_closure': 'PENDING_MC',
              'selection_version': VERSION, 'ROOT_version': ROOT.gROOT.GetVersion(), 'tables': []}
    for system in ('OO', 'pp'):
        dimension = 3 if system == 'OO' else 2
        names = [f'h{part}_{dimension}D' for part in ('Eff', 'Fak', 'Sec', 'Mul')]
        for selection in VERTEX:
            path = payload_path(system, selection)
            f = ROOT.TFile.Open(str(path), 'READ')
            if not f or f.IsZombie():
                raise RuntimeError('Cannot open tracking table: ' + str(path))
            try:
                hs = [f.Get(name) for name in names]
                if any(not h or not h.InheritsFrom('TH1') or h.GetDimension() != dimension for h in hs):
                    raise RuntimeError('Missing/incompatible histograms in ' + str(path))
                edges = axes(hs[0])
                if any(axes(h) != edges for h in hs[1:]):
                    raise RuntimeError('Tracking histogram axis mismatch: ' + str(path))
                if any(not all(math.isfinite(x) for x in edge) or any(b <= a for a, b in zip(edge, edge[1:])) for edge in edges):
                    raise RuntimeError('Invalid tracking axis edges')
                if edges[0][0] > -2.4 or edges[0][-1] < 2.4 or not edges[1][0] <= 2.0 < edges[1][-1]:
                    raise RuntimeError('Unexpected eta/pt axis coverage')
                if dimension == 3 and not edges[2][0] <= 0 <= 199 < edges[2][-1]:
                    raise RuntimeError('OO third axis does not cover hiBin 0..199')
                counts = {'regular_bins': 0, 'zero_or_low_efficiency': 0, 'nonfinite_content': 0,
                          'invalid_rate_content': 0, 'efficiency_above_one': 0}
                for ix in range(1, hs[0].GetNbinsX()+1):
                    for iy in range(1, hs[0].GetNbinsY()+1):
                        for iz in range(1, hs[0].GetNbinsZ()+1) if dimension == 3 else (0,):
                            index = hs[0].GetBin(ix, iy, iz)
                            efficiency, fake, secondary, multiple = [float(h.GetBinContent(index)) for h in hs]
                            counts['regular_bins'] += 1
                            if not all(math.isfinite(v) for v in (efficiency, fake, secondary, multiple)):
                                counts['nonfinite_content'] += 1
                            if efficiency <= .001: counts['zero_or_low_efficiency'] += 1
                            if efficiency > 1: counts['efficiency_above_one'] += 1
                            if fake < 0 or fake >= 1 or secondary < 0 or secondary >= 1 or multiple < 0:
                                counts['invalid_rate_content'] += 1
                report['tables'].append({'system': system, 'selection': selection, 'path': str(path.relative_to(REPO)),
                    'sha256': digest(path), 'bytes': path.stat().st_size, 'histograms': names,
                    'axes': edges, **counts, 'oo_lookup_pt_cap_GeV': 6.1 if system == 'OO' else None,
                    'note': 'Empty/invalid bins are counted, not filled or smoothed; event-level support is checked when MC is supplied.'})
            finally:
                f.Close()
    return report


def read_list(path):
    path = Path(path).resolve()
    result = []
    for line in path.read_text().splitlines():
        value = line.split('#',1)[0].strip()
        if not value: continue
        if '://' in value:
            raise ValueError('Stage XRootD inputs or use EOS-mounted paths before freezing this manifest: ' + value)
        item = Path(value).expanduser()
        item = (item if item.is_absolute() else path.parent/item).resolve()
        if not item.is_file(): raise FileNotFoundError(item)
        result.append(item)
    if not result or len(set(result)) != len(result):
        raise ValueError('Empty input list or repeated physical pathname: ' + str(path))
    return result


def inspect_input(path, system):
    import ROOT
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie(): raise RuntimeError('Cannot read input ' + str(path))
    try:
        event = f.Get('hiEvtAnalyzer/HiTree')
        if not event: raise ValueError('Missing hiEvtAnalyzer/HiTree')
        for name in ('run', 'lumi', 'evt', 'vz'):
            if not event.GetBranch(name): raise ValueError('Missing physical event field ' + name)
        tracks = f.Get('ppTrack/trackTree')
        if not tracks: tracks = f.Get('PbPbTracks/trackTree')
        if not tracks: raise ValueError('Missing track tree')
        for name in ('trkPt', 'trkEta', 'trkPhi', 'highPurity', 'trkPtError', 'trkCharge'):
            if not tracks.GetBranch(name): raise ValueError('Missing required track field ' + name)
        truth = f.Get('HiGenParticleAna/hi')
        if not truth: raise ValueError('Missing generator-particle tree HiGenParticleAna/hi')
        for name in ('pt', 'eta', 'phi', 'chg'):
            if not truth.GetBranch(name): raise ValueError('Missing generator field ' + name)
        jets = f.Get('akCs4PFJetAnalyzer/t' if system == 'OO' else 'ak4PFJetAnalyzer/t')
        if not jets: raise ValueError('Missing native R4 reconstructed jet tree')
        entries = int(event.GetEntries())
        if entries <= 0: raise ValueError('Empty input')
        if int(tracks.GetEntries()) != entries or int(truth.GetEntries()) != entries or int(jets.GetEntries()) != entries:
            raise ValueError('Event/track/truth/jet entry counts are not aligned')
        return {'path': str(path), 'entries': entries, 'bytes': path.stat().st_size,
                'root_uuid': str(f.GetUUID().AsString()), 'sha256': digest(path),
                'charge_class': tracks.GetBranch('trkCharge').GetClassName(),
                'truth_retention': 'Branch presence only; upstream generator-particle cuts still require the sample production information.'}
    finally:
        f.Close()


def plan(args):
    output = args.output.resolve()
    if output.exists(): raise FileExistsError('Use a fresh campaign directory: ' + str(output))
    samples = {}
    for system, filelist, namespace in (('OO', args.oo_list, args.oo_namespace), ('pp', args.pp_list, args.pp_namespace)):
        if filelist:
            if not namespace: raise ValueError('Provide the physical campaign namespace for ' + system)
            samples[system] = (namespace, [inspect_input(p, system) for p in read_list(filelist)])
    if not samples: raise ValueError('Supply at least one MC input list')
    report = audit_tables()
    output.mkdir(parents=True)
    cards = write_cards(output/'cards')
    (output/'tracking_tables.json').write_text(json.dumps(report, indent=2)+'\n')
    commands = ['#!/usr/bin/env bash', 'set -euo pipefail', 'cd ' + shlex.quote(str(REPO)), 'make -j2']
    def cmd(values): commands.append(shlex.join([str(v) for v in values]))
    python = sys.executable
    directory = REPO/'validation/run3Pairs'
    for selection in args.selections:
        folder = output/selection
        shard_dir = folder/'shards'; native_dir = folder/'native'
        shard_dir.mkdir(parents=True); native_dir.mkdir()
        tasks = []
        for system, (namespace, inputs) in samples.items():
            for index, item in enumerate(inputs):
                limit = item['entries'] if args.all_events else min(item['entries'], args.max_events)
                task = {'system': system, 'sample_namespace': namespace, 'file_index': index,
                        'input': item['path'], 'entry_limit': limit, 'input_identity': item,
                        'preparation_version': VERSION, 'track_selection': selection}
                tasks.append(task)
                prefix = f'{system}_{index:03d}'
                taskfile = folder/(prefix+'.task.json')
                taskfile.write_text(json.dumps(task, indent=2)+'\n')
                exported = native_dir/(prefix+'.jsonl')
                native = native_dir/(prefix+'.native')
                export = [python, directory/'export_forest.py', item['path'], exported,
                          '--system', system, '--sample-namespace', namespace, '--max-events', limit,
                          '--card', cards[(system,selection,'RecoReco')]]
                if index == 0: export += ['--parity', '--responses']
                cmd(export)
                cmd([python, directory/'batch_pairs.py', 'train', exported, shard_dir/prefix,
                     '--raw-ledger', native/'reco.jsonl.events.jsonl', '--metadata', str(exported)+'.metadata.json', '--task', taskfile])
        manifest = folder/'manifest.json'
        manifest.write_text(json.dumps({'version':VERSION,'selection':selection,'tasks':tasks,
                             'scope':'FULL_FILES' if args.all_events else 'BOUNDED_PREPARATION_PILOT'},indent=2)+'\n')
        for system in samples:
            training = folder/(system+'.training.json.gz'); result = folder/(system+'.closure.json.gz')
            cmd([python,directory/'batch_pairs.py','merge-training',manifest,shard_dir,training,
                 '--system',system,'--identity-db',folder/(system+'.identities.sqlite'),'--min-pairs',args.min_pairs])
            system_tasks = [t for t in tasks if t['system']==system]
            for task in system_tasks:
                prefix = f"{system}_{task['file_index']:03d}"
                cmd([python,directory/'batch_pairs.py','validate',shard_dir/(prefix+'.holdout.jsonl.gz'),training,shard_dir/(prefix+'.validation.json.gz')])
            cmd([python,directory/'batch_pairs.py','merge-validation',training,shard_dir,result])
            first = f"{system}_{system_tasks[0]['file_index']:03d}"
            cmd([python,directory/'batch_plot.py',result,folder/'plots'/system,'--first-then-all',
                 '--events',shard_dir/(first+'.holdout.jsonl.gz')])
    commandfile = output/'run_commands.sh'
    commandfile.write_text('\n'.join(commands)+'\n'); commandfile.chmod(0o755)
    print(json.dumps({'status':'PREPARED_NOT_EXECUTED','directory':str(output),'commands':str(commandfile),
                      'selections':args.selections,'systems':list(samples),'physics_closure':'PENDING_MC_EXECUTION'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    cards = sub.add_parser('cards'); cards.add_argument('--output',type=Path,required=True)
    audit = sub.add_parser('audit-tables'); audit.add_argument('--output',type=Path,required=True)
    p = sub.add_parser('plan'); p.add_argument('--output',type=Path,required=True)
    for system in ('oo','pp'):
        p.add_argument('--'+system+'-list',type=Path);p.add_argument('--'+system+'-namespace')
    p.add_argument('--selections',nargs='+',choices=list(VERTEX),default=['Nominal'])
    budget=p.add_mutually_exclusive_group();budget.add_argument('--all-events',action='store_true');budget.add_argument('--max-events',type=int,default=10000)
    p.add_argument('--min-pairs',type=int,default=10)
    args = parser.parse_args()
    if args.command == 'cards':
        result=write_cards(args.output);print(json.dumps({'cards':len(result),'output':str(args.output.resolve())}))
    elif args.command == 'audit-tables':
        result=audit_tables();args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    else:
        if args.max_events<1 or args.min_pairs<1: parser.error('Event/pair bounds must be positive')
        plan(args)


if __name__ == '__main__':
    main()
