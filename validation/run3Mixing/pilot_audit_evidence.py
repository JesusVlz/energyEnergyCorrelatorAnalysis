#!/usr/bin/env python3
"""Extra bounded-audit evidence: exact stored-particle equality and NH/vz support.

Coordinate support is before native reconstructed jet veto and jet selection;
it is an upper bound on final common-jet exposure, never a physics closure.
"""
import argparse
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def load(path):
    return [json.loads(line) for line in path.open()]


def coordinate_support(signals, donors, exclude_all_signal_backgrounds=False, exclude_all_signal_raw_ids=False):
    grid = defaultdict(list)
    all_signal_ids = {tuple(s['id']) for s in signals}
    all_signal_signatures = {s['background_identity'] for s in signals if s['background_identity']}
    for donor in donors:
        if not donor['pre_veto_selection_pass'] or donor['duplicate_physical_id']:
            continue
        if exclude_all_signal_backgrounds and donor['background_identity'] in all_signal_signatures:
            continue
        if exclude_all_signal_raw_ids and tuple(donor['id']) in all_signal_ids:
            continue
        grid[(donor['pool'], math.floor(donor['vz']/.5), donor['nh'])].append(donor)
    counts = Counter()
    by_activity = defaultdict(Counter)
    for signal in signals:
        if not signal['pre_veto_selection_pass'] or signal['duplicate_physical_id']:
            continue
        eligible = []
        tolerance = max(1., .01*signal['nh'])
        for pool in 'AB':
            found = False
            for dz, nh in itertools.product((-1,0,1), range(max(0,math.ceil(signal['nh']-tolerance)),math.floor(signal['nh']+tolerance)+1)):
                for donor in grid.get((pool, math.floor(signal['vz']/.5)+dz,nh), ()):
                    if donor['id'] == signal['id'] or abs(donor['vz']-signal['vz']) > .5:
                        continue
                    if signal['background_identity'] and signal['background_identity']==donor['background_identity']:
                        continue
                    found = True
                    break
                if found: break
            eligible.append(found)
        status = 'matched_A_and_B' if all(eligible) else 'no_'+('A' if not eligible[0] else '')+('B' if not eligible[1] else '')
        counts['signal_events_pre_veto'] += 1; counts[status] += 1
        hi = signal['hiBin']
        group = 'invalid' if not 0 <= hi < 200 else ('0_30' if hi < 60 else '30_50' if hi < 100 else '50_90' if hi < 180 else '90_100')
        by_activity[group]['signal_events_pre_veto'] += 1; by_activity[group][status] += 1
    total = counts['signal_events_pre_veto']
    return {'status': 'COORDINATE_FEASIBILITY_BEFORE_RECO_JET_VETO_AND_JET_SELECTION', 'counts':dict(counts),
            'fraction': counts['matched_A_and_B']/total if total else None,
            'activity': {key:dict(value) for key,value in sorted(by_activity.items())},
            'exclude_all_observed_signal_backgrounds_from_donors': exclude_all_signal_backgrounds,
            'exclude_all_scanned_signal_raw_ids_from_donors': exclude_all_signal_raw_ids,
            'limitations': 'No reconstructed veto or native signal jet selection here. Final analyzer coverage may be lower. No matching tolerance expansion.'}


def verify_exact_overlaps(matches):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    files = {}
    def particles(row, signal):
        path = row['source_file']
        if path not in files:
            source = ROOT.TFile.Open(path)
            tree = source.Get('HiGenParticleAna/hi')
            tree.SetBranchStatus('*',0)
            for name in ('pt','eta','phi','chg','sube'): tree.SetBranchStatus(name,1)
            files[path] = (source,tree)
        tree = files[path][1]
        tree.GetEntry(row['source_entry'])
        return sorted((float(pt),float(eta),float(phi),int(chg))
                      for pt,eta,phi,chg,sube in zip(tree.pt,tree.eta,tree.phi,tree.chg,tree.sube)
                      if not signal or sube > 0)
    counts = Counter()
    try:
        for match in matches:
            reference = particles(match['signals'][0],True)
            for signal in match['signals'][1:]:
                if particles(signal,True) != reference: raise ValueError('Fingerprint collision among signals')
            for donor in match['donors']:
                if particles(donor,False) != reference: raise ValueError('Fingerprint collision signal versus donor')
            counts['signatures_verified_exact_sorted_stored_particles'] += 1
            counts['signals_with_shared_stored_background'] += len(match['signals'])
            counts['donors_with_shared_stored_background'] += len(match['donors'])
            counts['different_raw_id_pairs'] += sum(s['id'] != d['id'] for s in match['signals'] for d in match['donors'])
    finally:
        for source,_ in files.values(): source.Close()
    return {'counts':dict(counts), 'comparison': 'exact sorted tuples of source stored float32 pt/eta/phi and integer charge; all signal sube>0 vs every pure-HIJING donor particle, including neutrals',
            'interpretation': 'Observed identical nonempty particle backgrounds across different physical event IDs; use observed background identity for exclusion/partition. Missing exact match cannot establish independent underlying lineage.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_dir',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists(): parser.error('Fresh output required')
    signals=load(args.input_dir/'signal.events.jsonl')
    donors=load(args.input_dir/'poolA.events.jsonl')+load(args.input_dir/'poolB.events.jsonl')
    result={'coordinate_support':coordinate_support(signals,donors),
            'coordinate_support_disjoint_from_all_observed_signal_backgrounds':coordinate_support(signals,donors,True),
            'coordinate_support_with_global_raw_ids_and_background_exclusion':coordinate_support(signals,donors,True,True),
            'shared_background_exact_verification':verify_exact_overlaps(load(args.input_dir/'background_identity_matches.jsonl'))}
    with args.output.open('x') as stream: json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
