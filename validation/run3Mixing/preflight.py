#!/usr/bin/env python3
"""Bounded measured-HF/vz pool feasibility study, never a background estimator.

Partner tolerances are explicit study settings. No calibration, activity weight,
jet-HF shift, component normalization or pair correction is inferred here.
"""
import argparse
import hashlib
import itertools
import json
import math
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'run3Pairs'))
from pair_counts import event_key, fold, physical_key
from export_forest import digest, select_branches


def match(events, vz_tolerance, hf_tolerance):
    """Deterministic two-partner assignment, isolated by physical-event fold.

    Match each signal event once, retaining the same partners for both radii and
    powers. Reuse is recorded, not treated as statistically independent exposure.
    """
    if not all(math.isfinite(x) and x > 0 for x in (vz_tolerance, hf_tolerance)):
        raise ValueError('Matching tolerances must be finite and positive')
    grid = defaultdict(list)
    seen = set()
    for event in events:
        key = physical_key(event)
        if (len(key) != 4 or not key[0] or key in seen or
                any(v < 0 or v >= 1 << 64 for v in key[1:])):
            raise ValueError('Invalid or duplicate physical event identity')
        seen.add(key)
        if any(not math.isfinite(event[k]) for k in ('vz', 'hf', 'generator_weight')):
            raise ValueError('Nonfinite event coordinate/weight')
        if event['hf'] < 0:
            raise ValueError('Negative measured HF')
        bucket = (key[0], fold(event), math.floor(event['vz']/vz_tolerance), math.floor(event['hf']/hf_tolerance))
        grid[bucket].append(event)
    assignments, failures, reuse = [], [], Counter()
    for signal in sorted(events, key=physical_key):
        if not signal['jets']:
            continue
        key = physical_key(signal)
        f, vzbin, hfbin = fold(signal), math.floor(signal['vz']/vz_tolerance), math.floor(signal['hf']/hf_tolerance)
        candidates = []
        for dv, dh in itertools.product((-1, 0, 1), repeat=2):
            for partner in grid.get((key[0], f, vzbin+dv, hfbin+dh), ()):
                if physical_key(partner) == key:
                    continue
                dz, dHF = partner['vz']-signal['vz'], partner['hf']-signal['hf']
                if abs(dz) <= vz_tolerance and abs(dHF) <= hf_tolerance:
                    candidates.append((abs(dz)/vz_tolerance+abs(dHF)/hf_tolerance,
                                       event_key(partner), physical_key(partner), partner))
        candidates.sort(key=lambda item: item[:3])
        if len(candidates) < 2:
            failures.append({'signal': list(key), 'fold': f, 'available_partners': len(candidates), 'jets': signal['jets']})
            continue
        partners = [item[3] for item in candidates[:2]]
        for partner in partners:
            reuse[physical_key(partner)] += 1
        assignments.append({'signal': list(key), 'fold': f, 'jets': signal['jets'],
                            'signal_vz': signal['vz'], 'signal_hf': signal['hf'],
                            'signal_hiBin': signal['activity'], 'signal_generator_weight': signal['generator_weight'],
                            'candidate_count': len(candidates),
                            'partners': [{'id': list(physical_key(p)), 'fold': fold(p), 'delta_vz': p['vz']-signal['vz'],
                                          'delta_hf': p['hf']-signal['hf'], 'hiBin': p['activity'],
                                          'generator_weight': p['generator_weight']} for p in partners]})
    coverage = {f'R{r:02d}_J{lo}_{hi}': {'signals': 0, 'matched': 0, 'failed': 0}
                for r in (4, 8) for lo, hi in ((40, 80), (80, 120), (120, 160))}
    for rows, state in ((assignments, 'matched'), (failures, 'failed')):
        for row in rows:
            for jet in row['jets']:
                for lo, hi in ((40, 80), (80, 120), (120, 160)):
                    if lo <= jet['pt'] < hi:
                        cell = coverage[f"R{int(round(jet['radius']*10)):02d}_J{lo}_{hi}"]
                        cell['signals'] += 1
                        cell[state] += 1
    hf_values = {e['hf'] for e in events}
    activity = {'min': min(hf_values, default=None), 'max': max(hf_values, default=None),
                'distinct_values': len(hf_values),
                'status': 'BLOCKED_DEGENERATE_COORDINATE' if len(hf_values) < 2 else 'VARIES_APPLICABILITY_NOT_VALIDATED'}
    return {'status': 'DIAGNOSTIC_ONLY_NOT_VALIDATED', 'activity_coordinate': activity, 'vz_tolerance_cm': vz_tolerance,
            'hf_tolerance_GeV': hf_tolerance, 'partner_rule': 'two distinct other physical events, same sample namespace and fold, nearest normalized L1 mismatch; FNV identity tie-break',
            'pool_events': len(events), 'signal_events': len(assignments)+len(failures),
            'matched_signal_events': len(assignments), 'failed_signal_events': len(failures),
            'coverage_jets_shared_by_both_powers': coverage, 'assignments': assignments, 'failures': failures,
            'partner_reuse': [{'id': list(k), 'uses': v} for k, v in sorted(reuse.items())],
            'max_partner_reuse': max(reuse.values(), default=0),
            'covariance': 'NOT_RUN; signal/partner identities and signed weights retained; reuse must propagate jointly'}


def scan(args):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    file = ROOT.TFile.Open(str(args.input))
    if not file or file.IsZombie():
        raise ValueError('Cannot open forest')
    event, skim = file.Get('hiEvtAnalyzer/HiTree'), file.Get('skimanalysis/HltTree')
    if not event or not skim:
        raise ValueError('Missing event/skim tree')
    branches = ['run', 'lumi', 'evt', 'vz', 'hiHF_pf', 'hiBin', 'weight']
    select_branches(event, branches)
    filters = (['pprimaryVertexFilter', 'OOphfCoincFilterPF2Th4', 'pileupVertexFilter']
               if args.system == 'OO' else ['pprimaryVertexFilter', 'ppNoScrapingFilter'])
    select_branches(skim, filters)
    native, missing = {}, {}
    for radius in (.4, .8):
        name = ('akCs' if args.system == 'OO' else 'ak')+str(round(radius*10))+'PFJetAnalyzer/t'
        tree = file.Get(name)
        if not tree:
            missing[str(radius)] = name
            continue
        select_branches(tree, ['nref', 'jtpt', 'rawpt', 'WTAeta', 'trackMax'])
        native[radius] = tree
    nentries = int(event.GetEntries())
    if any(int(t.GetEntries()) != nentries for t in (skim, *native.values())):
        raise ValueError('Unaligned forest trees')
    counts = Counter({key: 0 for key in ('scanned', 'selection_rejected',
                     'source_hiBin_invalid', 'source_hiBin_0_90', 'source_hiBin_90_100',
                     'accepted_hiBin_invalid', 'accepted_hiBin_0_90', 'accepted_hiBin_90_100')})
    events = []
    seen = set()
    for entry in range(min(args.events, nentries)):
        event.GetEntry(entry)
        skim.GetEntry(entry)
        counts['scanned'] += 1
        flags = {name: bool(getattr(skim, name)) for name in filters}
        for name, passed in flags.items():
            counts['passes_'+name] += passed
        hibin = int(event.hiBin)
        group = 'invalid' if not 0 <= hibin < 200 else ('90_100' if hibin >= 180 else '0_90')
        counts['source_hiBin_'+group] += 1
        if abs(float(event.vz)) > 15 or not all(flags.values()):
            counts['selection_rejected'] += 1
            continue
        counts['accepted_hiBin_'+group] += 1
        row = {'sample': args.system, 'sample_namespace': args.sample_namespace,
               'id': [int(event.run), int(event.lumi), int(event.evt)], 'source_entry': entry,
               'vz': float(event.vz), 'hf': float(event.hiHF_pf), 'activity': hibin,
               'generator_weight': float(event.weight), 'jets': []}
        key = physical_key(row)
        if key in seen:
            raise ValueError('Duplicate physical ID in selected forest')
        seen.add(key)
        for radius, tree in native.items():
            tree.GetEntry(entry)
            for j in range(int(tree.nref)):
                raw = float(tree.rawpt[j])
                if 40 <= float(tree.jtpt[j]) < 160 and abs(float(tree.WTAeta[j])) < 1.6 and raw > 0 and .01 < float(tree.trackMax[j])/raw < .98:
                    row['jets'].append({'radius': radius, 'index': j, 'pt': float(tree.jtpt[j])})
        events.append(row)
    file.Close()
    return events, {'counts': dict(counts), 'native_tree_missing': missing, 'total_source_entries': nentries,
                    'entry_range': [0, min(args.events, nentries)], 'sample_namespace': args.sample_namespace,
                    'selection': 'abs(vz)<=15 && '+' && '.join(filters)+'; no MC trigger; no hiBin cut',
                    'input_sha256': digest(args.input), 'input': str(args.input.resolve()),
                    'jet_selection': 'native stored jtpt in [40,160),abs(WTAeta)<1.6,.01<trackMax/rawpt<.98; JEC applicability unresolved',
                    'source_bias': 'local single-file MC only; upstream skim loss, MB representativeness and data/MC HF compatibility unmeasured',
                    'root_version': ROOT.gROOT.GetVersion()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--system', choices=('OO', 'pp'), required=True)
    parser.add_argument('--sample-namespace', required=True)
    parser.add_argument('--events', type=int, default=1000)
    parser.add_argument('--vz-tolerance', type=float, required=True)
    parser.add_argument('--hf-tolerance', type=float, required=True)
    args = parser.parse_args()
    if not 1 <= args.events <= 10000:
        parser.error('Bounded scan requires 1..10000 events')
    start = time.monotonic()
    events, metadata = scan(args)
    result = match(events, args.vz_tolerance, args.hf_tolerance)
    ledger = args.output.with_suffix('.events.jsonl')
    with ledger.open('x') as stream:
        for event in events:
            stream.write(json.dumps(event, allow_nan=False)+'\n')
    result.update(metadata)
    result.update({'ledger': str(ledger.resolve()), 'ledger_sha256': digest(ledger),
                   'tool_sha256': digest(Path(__file__)), 'identity_tool_sha256': digest(Path(__file__).resolve().parents[1]/'run3Pairs/pair_counts.py'),
                   'elapsed_seconds': time.monotonic()-start,
                   'physics_blockers': ['HF/vz tolerances are study settings, not validated matching prescription',
                                        'jet contribution to HF and matching bias NOT_VALIDATED',
                                        'no independent data/MB partner source or background component closure',
                                        'source-filter and calibration applicability unresolved']})
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('status', 'pool_events', 'signal_events', 'matched_signal_events', 'failed_signal_events', 'max_partner_reuse')}))


if __name__ == '__main__':
    main()
