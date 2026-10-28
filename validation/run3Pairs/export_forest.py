#!/usr/bin/env python3
"""Bounded real-file diagnostic exporter; never manufactures native R8 inputs.

Uses native forest WTA axes and stored jtpt (not rederived JEC), mirrors the
current Run-3 nominal track cuts, and reproduces inherited relative-eta tracking
lookups. These explicit limitations BLOCK physics promotion. Original detector
coordinates and physical-event identities are retained. No analyzer is altered.
"""
import argparse
import hashlib
import json
import math
import time
from collections import Counter
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def charged(value):
    return (ord(value) if isinstance(value, str) else int(value)) != 0


def select_branches(tree, branches):
    missing = [name for name in branches if not tree.GetBranch(name)]
    if missing:
        raise ValueError(str(tree.GetName()) + ' missing required branches: ' + ','.join(missing))
    tree.SetBranchStatus('*', 0)
    for name in branches:
        tree.SetBranchStatus(name, 1)


def selected_reco(tree, i):
    # Mirrors EECAnalyzer::PassTrackCuts Run-3 branch, strict pT>2 for this
    # measurement. No unsupported legacy calo/chi2/hit cuts are introduced.
    pt, eta = float(tree.trkPt[i]), float(tree.trkEta[i])
    dzerr, dxyerr = float(tree.trkDzErrFirstVtx[i]), float(tree.trkDxyErrFirstVtx[i])
    return (charged(tree.trkCharge[i]) and 2 < pt < 500 and abs(eta) < 2.4 and
            bool(tree.highPurity[i]) and float(tree.trkPtError[i])/pt < .1 and
            dzerr > 0 and dxyerr > 0 and abs(float(tree.trkDzFirstVtx[i])/dzerr) < 3 and
            abs(float(tree.trkDxyFirstVtx[i])/dxyerr) < 3)


def lookup(payload, system, pt, eta, hibin):
    if system == 'OO' and not 0 <= hibin <= 199:
        raise ValueError('Invalid OO hiBin; no edge clamping')
    coordinates = (eta, min(pt, 6.1), hibin) if system == 'OO' else (eta, pt)
    values = []
    for h in payload:
        axes = [h.GetXaxis(), h.GetYaxis(), h.GetZaxis()][:h.GetDimension()]
        if any(not math.isfinite(v) or v < a.GetXmin() or v >= a.GetXmax() for v, a in zip(coordinates, axes)):
            raise ValueError('Out-of-domain single-track lookup; no overflow/unity fallback')
        values.append(float(h.GetBinContent(h.FindBin(*coordinates))))
    eff, fake, multiple, secondary = values
    if not all(math.isfinite(x) for x in values) or eff <= .001 or eff*(1+multiple) <= 0:
        raise ValueError('Invalid single-track payload content')
    correction = (1-fake)*(1-secondary)/(eff*(1+multiple))
    if not math.isfinite(correction) or correction <= 0:
        raise ValueError('Invalid single-track correction')
    return correction


def export(args):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    started = time.monotonic()
    f = ROOT.TFile.Open(str(args.input))
    if not f or f.IsZombie():
        raise ValueError('Cannot open input forest')
    event = f.Get('hiEvtAnalyzer/HiTree')
    truth = f.Get('HiGenParticleAna/hi')
    tracks = f.Get('PbPbTracks/trackTree' if args.system == 'OO' else 'ppTracks/trackTree')
    if not event or not truth or not tracks:
        raise ValueError('Missing event/truth/track tree')
    select_branches(event, ['run', 'lumi', 'evt', 'vz', 'weight', 'hiBin'])
    select_branches(truth, ['pt', 'eta', 'phi', 'chg', 'sube'])
    select_branches(tracks, ['trkPt', 'trkEta', 'trkPhi', 'trkCharge', 'highPurity', 'trkPtError',
                             'trkDzFirstVtx', 'trkDzErrFirstVtx', 'trkDxyFirstVtx', 'trkDxyErrFirstVtx'])
    counts = Counter({name: 0 for name in ('scanned_events', 'exported_events', 'rejected_vz', 'unsupported_hiBin_events', 'selected_jets_R04', 'selected_jets_R08')})
    payload_file = ROOT.TFile.Open(str(args.payload))
    if not payload_file or payload_file.IsZombie():
        raise ValueError('Cannot open nominal single-track payload')
    names = ['hEff_3D', 'hFak_3D', 'hMul_3D', 'hSec_3D'] if args.system == 'OO' else ['hEff_2D', 'hFak_2D', 'hMul_2D', 'hSec_2D']
    payload = [payload_file.Get(name) for name in names]
    if not all(payload):
        raise ValueError('Missing nominal single-track histogram')
    native = {}
    blocked = {}
    radius_seconds = {}
    for radius in (.4, .8):
        name = ('akCs' if args.system == 'OO' else 'ak') + ('4' if radius == .4 else '8') + 'PFJetAnalyzer/t'
        tree = f.Get(name)
        if not tree:
            blocked[str(radius)] = 'BLOCKED_NATIVE_TREE_MISSING: ' + name
            continue
        select_branches(tree, ['nref', 'jtpt', 'rawpt', 'WTAeta', 'WTAphi', 'trackMax', 'refpt',
                               'ngen', 'genpt', 'WTAgeneta', 'WTAgenphi', 'genmatchindex'])
        native[radius] = tree
        radius_seconds[str(radius)] = 0.
    n = min(args.max_events, int(event.GetEntries()))
    expected_entries = int(event.GetEntries())
    if any(int(t.GetEntries()) != expected_entries for t in [truth, tracks, *native.values()]):
        raise ValueError('Unaligned event/track/truth/native-jet entry counts')
    identities = set()
    with args.output.open('w') as output:
        for entry in range(n):
            event.GetEntry(entry)
            counts['scanned_events'] += 1
            if abs(float(event.vz)) > 15:
                counts['rejected_vz'] += 1
                continue
            hibin = int(event.hiBin)
            if args.system == 'OO' and not 0 <= hibin <= 199:
                counts['unsupported_hiBin_events'] += 1
                continue
            truth.GetEntry(entry)
            tracks.GetEntry(entry)
            identity = (int(event.run), int(event.lumi), int(event.evt))
            if identity in identities:
                raise ValueError('Duplicate physical event ID within input')
            identities.add(identity)
            row = {'sample': args.system, 'sample_namespace': args.sample_namespace,
                   'id': list(identity), 'source_entry': entry, 'activity': hibin,
                   'generator_weight': float(event.weight), 'event_weight': 1., 'jets': []}
            selected = [i for i in range(len(tracks.trkPt)) if selected_reco(tracks, i)]
            generated = [i for i in range(len(truth.pt)) if charged(truth.chg[i]) and 2 < float(truth.pt[i]) < 500 and abs(float(truth.eta[i])) < 2.4]
            for radius, jets in native.items():
                radius_started = time.monotonic()
                jets.GetEntry(entry)
                for j in range(int(jets.nref)):
                    pt, eta, phi, raw = float(jets.jtpt[j]), float(jets.WTAeta[j]), float(jets.WTAphi[j]), float(jets.rawpt[j])
                    if not (40 <= pt < 160 and abs(eta) < 1.6 and raw > 0 and .01 < float(jets.trackMax[j])/raw < .98):
                        continue
                    jet = {'radius': radius, 'index': j, 'pt': pt, 'eta': eta, 'phi': phi,
                           'refpt': float(jets.refpt[j]), 'reco': [], 'gen': []}
                    # Reco and gen pair references use exactly this selected reco
                    # jet, without switching to a matched generator axis/population.
                    for population, tree, indices in [('reco', tracks, selected), ('gen', truth, generated)]:
                        for index in indices:
                            prefix = 'trk' if population == 'reco' else ''
                            tpt = float(getattr(tree, prefix+'Pt' if prefix else 'pt')[index])
                            teta = float(getattr(tree, prefix+'Eta' if prefix else 'eta')[index])
                            tphi = float(getattr(tree, prefix+'Phi' if prefix else 'phi')[index])
                            if math.hypot(teta-eta, math.remainder(tphi-phi, 2*math.pi)) >= radius:
                                continue
                            track = {'index': index, 'pt': tpt, 'eta': teta, 'phi': tphi}
                            if population == 'reco':
                                track['correction'] = lookup(payload, args.system, tpt, teta-eta, hibin)
                                track['single_track_lookup_eta'] = teta-eta
                            else:
                                track['subevent'] = int(truth.sube[index])
                            jet[population].append(track)
                    row['jets'].append(jet)
                    counts['selected_jets_R04' if radius == .4 else 'selected_jets_R08'] += 1
                radius_seconds[str(radius)] += time.monotonic() - radius_started
            output.write(json.dumps(row, allow_nan=False) + '\n')
            counts['exported_events'] += 1
    metadata = {'status': 'BLOCKED_DIAGNOSTIC_ONLY', 'input': str(args.input.resolve()), 'input_sha256': digest(args.input),
                'payload': str(args.payload.resolve()), 'payload_sha256': digest(args.payload), 'output_sha256': digest(args.output),
                'system': args.system, 'sample_namespace': args.sample_namespace, 'entry_range': [0, n],
                'counts': dict(counts), 'native_radius_blockers': blocked,
                'radius_jet_track_loop_seconds': radius_seconds,
                'timing_scope': 'per-radius jet/track loop; shared I/O, selection, hashing and serialization included only in total elapsed time',
                'selection': 'vz<=15; strict Run-3 nominal track cuts, pt>2,abs(eta)<2.4; WTA abs(eta)<1.6; stored forest jtpt 40..160; trackMax/rawpt .01..98',
                'source_identity': 'sample_namespace + run,lumi,evt; event and track source index retained',
                'truth_selection': 'all stored charged final-particle collection entries, pT>2,abs(eta)<2.4; all subevents (SubeventCut2); final status follows inherited GeneratorLevelForestReader::GetTrackMCStatus=1',
                'input_threshold_requirement': 'separate low-pT audit required; no 2GeV filename assumption',
                'single_track_coordinates': 'inherited relative eta=teta-WTAjeteta reproduced; original detector eta/phi retained; no convention change',
                'event_weight_policy': 'stored generator weight once; extra event weight=1; no unsupported centrality/vz replacement or extrapolation',
                'blockers': ['stored forest jtpt: current applicable JEC not established by this exporter', 'full MC skim/event filter equivalence to analyzer NOT_VALIDATED',
                             'single-track calibration/reconstruction applicability not established from filename or axes',
                             *(['OO hiBin validity awaits experts'] if args.system == 'OO' else []),
                             'large-angle control and n1/n2 weighted closure NOT_VALIDATED'],
                'tool_sha256': digest(Path(__file__)),
                'track_cut_source_sha256': digest(Path(__file__).resolve().parents[2] / 'src/EECAnalyzer.cxx'),
                'elapsed_seconds': time.monotonic()-started, 'root_version': ROOT.gROOT.GetVersion()}
    args.output.with_suffix(args.output.suffix+'.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))
    payload_file.Close()
    f.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--system', choices=('OO', 'pp'), required=True)
    parser.add_argument('--sample-namespace', required=True, help='same physical campaign across radius streams, distinct for independent samples')
    parser.add_argument('--payload', type=Path, required=True)
    parser.add_argument('--max-events', type=int, default=100)
    args = parser.parse_args()
    if not 1 <= args.max_events <= 10000:
        parser.error('bounded exporter requires 1..10000 events')
    export(args)


if __name__ == '__main__':
    main()
