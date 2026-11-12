#!/usr/bin/env python3
"""Export R4 Run-3 pairs with the analyzer's selections, JEC and corrections.

Two native RecoReco/RecoGen runs on identical entries supply a common reco-jet
population. Their events, jets and weights must match exactly before merging.
The retained native ROOT files permit a bounded JSON/EEC sum comparison.
"""
import argparse
import hashlib
import itertools
import json
import math
import subprocess
import time
from array import array
from pathlib import Path
from preparation_contract import selection_from_card, payload_path, export_contract

REPO = Path(__file__).resolve().parents[2]
PAYLOADS = {
    'OO': REPO/'trackCorrectionTables/OO2025/Eff_OO_2025_Hijing_MB_Centrality_fromHihfpf_NoPU_3D_Nominal_Official_18Nov2025.root',
    'pp': REPO/'trackCorrectionTables/pp2024/Eff_ppref_2024_Pythia_QCD_pThat15_NopU_2D_Nominal_Official.root',
}


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def card_text(base, settings):
    lines = []
    text = base.read_text()
    baseline = {}
    for line in text.splitlines():
        words = line.split('#', 1)[0].split()
        if len(words) > 1: baseline[words[0]] = words[1]
    run3 = int(float(settings.get('DataType', baseline.get('DataType', -1)))) in (9, 10, 11, 12)
    unused = {'CalorimeterSignalLimitPt', 'HighPtEtFraction', 'Chi2QualityCut', 'MinimumTrackHits'} if run3 else set()
    remaining = {key: value for key, value in settings.items() if key not in unused}
    for line in text.splitlines():
        fields = line.split()
        key = fields[0] if fields else ''
        if key in unused: continue
        if key in settings:
            if key in remaining:
                lines.append(f'{key} {remaining.pop(key)}')
        else:
            lines.append(line)
    lines.extend(f'{key} {value}' for key, value in remaining.items())
    return '\n'.join(lines)+'\n'


def rows(path):
    with path.open() as stream:
        for line in stream:
            yield json.loads(line)


def parity(root_path, expected):
    import ROOT
    f = ROOT.TFile.Open(str(root_path))
    sparse = f.Get('energyEnergyCorrelator')
    if not sparse:
        raise ValueError('Native energyEnergyCorrelator histogram missing')
    observed = [0., 0.]
    coordinates = array('i', [0]*sparse.GetNdimensions())
    for index in range(sparse.GetNbins()):
        value = sparse.GetBinContent(index, coordinates)
        if coordinates[4] != 1:
            continue
        exponent = coordinates[6]-1
        if exponent in (0, 1):
            observed[exponent] += value
    checks = []
    for exponent in (1, 2):
        want = expected[exponent-1]
        got = observed[exponent-1]
        relative = abs(got-want)/max(abs(want), 1e-30)
        checks.append({'n': exponent, 'native': got, 'exported': want, 'relative_difference': relative,
                       'pass': relative < 2e-5 or abs(got-want) < 1e-10})
    if not all(check['pass'] for check in checks):
        raise ValueError('Native/exported EEC weight mismatch: '+str(checks))
    f.Close()
    return checks


def export(args):
    started = time.monotonic()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent/(output.stem+'.native')
    work.mkdir(exist_ok=True)
    base = args.card or REPO/('cardEECOOMCInclusiveR04.input' if args.system == 'OO' else 'cardEECpp5v36TeVMCR04.input')
    selection = selection_from_card(base)
    payload = payload_path(args.system, selection)
    contract = export_contract(args.system, selection, payload)
    if args.payload and digest(args.payload) != digest(payload):
        raise ValueError('--payload differs from the analyzer nominal payload; configure analyzer explicitly')
    base = args.card or REPO/('cardEECOOMCInclusiveR04.input' if args.system == 'OO' else 'cardEECpp5v36TeVMCR04.input')
    if args.reuse_native:
        raise ValueError('Native reuse requires a verified campaign receipt; regenerate this standalone export after the Run-3 cut repair')
    paths, commands = {}, []
    for population, correlation in [('reco', 0), ('gen', 1)]:
        prefix = work/population
        paths[population] = prefix.with_suffix('.jsonl')
        settings = {'McCorrelationType': correlation, 'BackgroundMethods': 0,
                    'MinJetPtCut': 40, 'MaxJetPtCut': 160, 'MinTrackPtCut': 0.5,
                    'TrackEtaCut': 2.4, 'WeightExponent': '1 2', 'MinPairTrackPt': 2,
                    'MaxEventsPerFile': args.max_events, 'HFEnergyShift': 0,
                    'FilledHistograms': 279 if args.responses and population == 'reco' else 23,
                    'DisableTrackPairEfficiencyCorrection': 1, 'JetPtWeight': 0,
                    'SkipCovarianceMatrix': 1, 'OOClosurePilot': 0,
                    'CentralityBinEdges': '0 30 50 90 100' if args.system == 'OO' else '-1 100',
                    'JetPtBinEdgesEEC': '40 80 120 160',
                    'PairExportFile': paths[population], 'PairExportNamespace': args.sample_namespace}
        card = prefix.with_suffix('.input')
        card.write_text(card_text(base, settings))
        command = [str(REPO/'eecAnalysis'), str(args.input.resolve()), str(card),
                   str(prefix.with_suffix('.root')), '0', 'true']
        commands.append(command)
        if not args.reuse_native:
            with prefix.with_suffix('.log').open('w') as log:
                subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, check=True)
    if digest(Path(str(paths['reco'])+'.events.jsonl')) != digest(Path(str(paths['gen'])+'.events.jsonl')):
        raise ValueError('Native RecoReco/RecoGen scanned-event ledgers differ')
    counts = {'exported_events': 0, 'selected_jets_R04': 0, 'reco_tracks': 0, 'gen_tracks': 0,
              'factor_changed_from_relative_eta': 0, 'jet_pt_changed_from_stored': 0,
              'unsupported_auxiliary_reco_tracks': 0}
    expected = {'reco': [0., 0.], 'gen': [0., 0.]}
    max_factor_change = 0.
    identities = set()
    with output.open('w') as result:
        for reco, gen in itertools.zip_longest(rows(paths['reco']), rows(paths['gen'])):
            if reco is None or gen is None:
                raise ValueError('Native selected-event counts differ')
            for key in ('sample_namespace', 'source_file', 'source_entry', 'id', 'activity', 'generator_weight', 'event_weight'):
                if reco[key] != gen[key]:
                    raise ValueError('Native event selection/weight mismatch: '+key)
            identity = (reco['sample_namespace'], *reco['id'])
            if identity in identities:
                raise ValueError('Duplicate namespaced physical event identity')
            identities.add(identity)
            reco['reco'] = reco.pop('tracks')
            reco['gen'] = gen.pop('tracks')
            counts['exported_events'] += 1
            counts['reco_tracks'] += len(reco['reco'])
            counts['gen_tracks'] += len(reco['gen'])
            counts['unsupported_auxiliary_reco_tracks'] += sum(bool(track.get('unsupported_single')) for track in reco['reco'])
            for rjet, gjet in itertools.zip_longest(reco['jets'], gen['jets']):
                if rjet is None or gjet is None:
                    raise ValueError('Native common-reco-jet counts differ')
                for key in ('radius', 'index', 'pt', 'eta', 'phi', 'jet_weight'):
                    if rjet[key] != gjet[key]:
                        raise ValueError('Native common-reco-jet selection mismatch: '+key)
                rjet['reco'] = rjet.pop('tracks')
                rjet['gen'] = gjet.pop('tracks')
                counts['selected_jets_R04'] += 1
                counts['jet_pt_changed_from_stored'] += not math.isclose(rjet['pt'], rjet['stored_pt'], rel_tol=1e-6)
                for track in rjet['reco']:
                    old_factor = track['legacy_relative_eta_correction']
                    relative = abs(track['correction']/old_factor-1) if math.isfinite(old_factor) and old_factor > 0 else 1.
                    counts['factor_changed_from_relative_eta'] += relative > 1e-6
                    max_factor_change = max(max_factor_change, relative)
                if args.parity:
                    for population in ('reco', 'gen'):
                        for first, second in itertools.combinations(rjet[population], 2):
                            weight = reco['event_weight']*reco['generator_weight']*rjet['jet_weight']
                            weight *= first['correction']*second['correction']
                            momentum = first['pt']*second['pt']
                            for exponent in (1, 2):
                                expected[population][exponent-1] += weight*momentum**exponent
            reco['tracking_contract'] = contract
            result.write(json.dumps(reco, allow_nan=False)+'\n')
    native_parity = {population: parity(work/(population+'.root'), expected[population])
                     for population in ('reco', 'gen')} if args.parity else 'not requested (reuse bounded parity evidence)'
    metadata = {'status': 'NATIVE_SELECTION_EXPORT', 'system': args.system, 'tracking_contract': contract,
                'sample_namespace': args.sample_namespace, 'input': str(args.input.resolve()),
                'input_size_bytes': args.input.stat().st_size,
                'payload': str(payload), 'payload_sha256': digest(payload),
                'entry_limit_per_file': args.max_events, 'counts': counts,
                'commands': commands, 'base_card': str(base), 'native_outputs': str(work),
                'native_selection_weight_identity': 'PASS: exact selected events, common reconstructed jets and event/jet weights in RecoReco and RecoGen',
                'native_eec_parity': native_parity,
                'maximum_single_factor_change_relative_eta': max_factor_change,
                'single_track_coordinates': 'physical detector eta and source-event hiBin; OO pT capped at 6.1 GeV by native corrector',
                'event_weight_policy': 'native Run-3 generator factor=1; stored forest weight retained separately; native vz and OO centrality fits',
                'jet_policy': 'native Run-3 R4 raw pT + applicable L2Relative; nominal card Jet ID and whole-event veto retained',
                'unsupported': ['no cross-section generator weight calibration inferred; analyzer existing unity policy retained',
                                'auxiliary full-event tracks with unavailable single-track factors are flagged; any consumed all-event pair cell/control must be invalidated'],
                'tool_sha256': digest(Path(__file__)), 'analyzer_sha256': digest(REPO/'src/EECAnalyzer.cxx'), 'native_binary_sha256': digest(REPO/'eecAnalysis'),
                'output_sha256': digest(output), 'elapsed_seconds': time.monotonic()-started}
    output.with_suffix(output.suffix+'.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--system', choices=('OO', 'pp'), required=True)
    parser.add_argument('--sample-namespace', required=True)
    parser.add_argument('--payload', type=Path)
    parser.add_argument('--max-events', type=int, default=10000)
    parser.add_argument('--card', type=Path)
    parser.add_argument('--responses', action='store_true', help='also fill existing native AN angular/momentum response histograms')
    parser.add_argument('--parity', action='store_true', help='bounded native/exported EEC weight comparison, both powers')
    parser.add_argument('--reuse-native', action='store_true', help='merge/check already completed native runs')
    args = parser.parse_args()
    if args.max_events < 1:
        parser.error('--max-events must be positive')
    export(args)


if __name__ == '__main__':
    main()
