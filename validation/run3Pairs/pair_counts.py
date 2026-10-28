#!/usr/bin/env python3
"""Candidate residual pair maps and disjoint-event weighted closure.
No output from this module is a production payload. Energy weighting is excluded
from map derivation. A map multiplies each n=1/n=2 pair exactly once.
"""
import argparse
import bisect
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

PT = (2., 3., 4., 8., 12., 20., 50., 100., 500.)
JET = (40., 80., 120., 160.)
ACTIVITY = (0., 60., 100., 180., 200.)  # stored hiBin, provisional activity only
DR = (0., .002, .005, .01, .02, .04, .07, .1, .15, .2, .3, .4, .6, .8, 1., 1.2, 1.6)
CONTROL = (.1, .4)  # inherited starting control, NOT a validated physics range


def physical_key(event):
    """No radius/power/file/entry in physical identity; samples stay distinct."""
    return (str(event.get('sample_namespace', event['sample'])), *(int(v) for v in event['id']))


def event_key(event):
    # Byte-identical to src/OOEventIdentity.h::EventKey; physical namespace is
    # required because independent samples can reuse numerical run/lumi/evt.
    key = physical_key(event)
    value = 14695981039346656037
    data = key[0].encode('utf-8') + b'\0' + b''.join(v.to_bytes(8, 'little', signed=False) for v in key[1:])
    for byte in data:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return value


def fold(event):
    return event_key(event) % 2


def bin_index(value, edges):
    if not math.isfinite(value) or value < edges[0] or value >= edges[-1]:
        return None
    return bisect.bisect_right(edges, value) - 1


def separation(a, b):
    return math.hypot(a['eta'] - b['eta'], math.remainder(a['phi'] - b['phi'], 2 * math.pi))


def pair_weight(a, b, n, event_weight, correction=1.):
    return event_weight * a.get('correction', 1.) * b.get('correction', 1.) * (a['pt'] * b['pt']) ** n * correction


def pairs(tracks):
    seen = set()
    for track in tracks:
        if track['index'] in seen:
            raise ValueError('Duplicate physical track index within one population')
        seen.add(track['index'])
    return itertools.combinations(tracks, 2)


def base_key(event, jet, a, b):
    pt_hi, pt_lo = sorted((a['pt'], b['pt']), reverse=True)
    activity = bin_index(event['activity'], ACTIVITY) if event['sample'] == 'OO' else 0
    indices = (activity, bin_index(jet['pt'], JET), bin_index(pt_hi, PT), bin_index(pt_lo, PT))
    if any(i is None for i in indices):
        return None
    return (event['sample'], 'R04' if jet['radius'] == .4 else 'R08', *indices)


def validate(events):
    seen = set()
    for event in events:
        identity = physical_key(event)
        if len(identity) != 4 or identity in seen:
            raise ValueError('Missing/duplicate physical event identity')
        seen.add(identity)
        if event['sample'] not in ('OO', 'pp') or not identity[0] or any(v < 0 or v >= 1 << 64 for v in identity[1:]):
            raise ValueError('Unsupported system or invalid physical event identity')
        if event['sample'] == 'OO' and bin_index(event['activity'], ACTIVITY) is None:
            raise ValueError('Unsupported OO activity; no bin clamping or silent exclusion')
        for key in ('generator_weight', 'event_weight'):
            if not math.isfinite(event[key]):
                raise ValueError('Nonfinite event weight')
        for jet in event['jets']:
            if bin_index(jet['pt'], JET) is None:
                raise ValueError('Unsupported reconstructed jet pT')
            if jet['radius'] not in (.4, .8):
                raise ValueError('Unsupported native radius')
            for population in ('reco', 'gen'):
                for track in jet[population]:
                    if not (2 < track['pt'] < 500 and abs(track['eta']) < 2.4 and math.isfinite(track['phi'])):
                        raise ValueError('Input violates strict measurement track acceptance')
                    if not math.isfinite(track.get('correction', 1.)) or track.get('correction', 1.) <= 0:
                        raise ValueError('Invalid single-track correction; no unity fallback')
                list(pairs(jet[population]))


def derive(events, min_pairs=10, control=CONTROL):
    """n=0 pair counts, with generator/event and nominal single-track factors.

    Sumw2 is bookkeeping only. Correlated pair/event errors require replicas;
    per-event closure vectors are retained and no independent-Poisson claim made.
    """
    cells = defaultdict(lambda: {'gen': 0., 'single': 0., 'raw': 0., 'gen_sumw2': 0., 'single_sumw2': 0., 'gen_pairs': 0, 'reco_pairs': 0})
    for event in events:
        if fold(event) != 0:
            continue
        weight = event['generator_weight'] * event['event_weight']
        for jet in event['jets']:
            for population in ('gen', 'reco'):
                for a, b in pairs(jet[population]):
                    base = base_key(event, jet, a, b)
                    dr = bin_index(separation(a, b), DR)
                    if base is None or dr is None:
                        continue
                    cell = cells[(*base, dr)]
                    if population == 'gen':
                        cell['gen'] += weight
                        cell['gen_sumw2'] += weight * weight
                        cell['gen_pairs'] += 1
                    else:
                        corrected = pair_weight(a, b, 0, weight)
                        cell['raw'] += weight
                        cell['single'] += corrected
                        cell['single_sumw2'] += corrected * corrected
                        cell['reco_pairs'] += 1
    controls = defaultdict(lambda: {'gen': 0., 'single': 0., 'gen_pairs': 0, 'reco_pairs': 0})
    for key, cell in cells.items():
        dr = key[-1]
        if DR[dr] >= control[0] and DR[dr+1] <= control[1]:
            for name in controls[key[:-1]]:
                controls[key[:-1]][name] += cell[name]
    result = {}
    for key, cell in cells.items():
        large = controls[key[:-1]]
        alpha = large['single'] / large['gen'] if large['gen'] > 0 else None
        valid = (cell['gen_pairs'] >= min_pairs and cell['reco_pairs'] >= min_pairs and
                 large['gen_pairs'] >= min_pairs and large['reco_pairs'] >= min_pairs and
                 cell['gen'] > 0 and cell['single'] > 0 and alpha is not None and alpha > 0)
        correction = alpha * cell['gen'] / cell['single'] if valid else None
        result[key] = dict(cell, control=dict(large), large_angle_single_over_gen=alpha,
                           correction=correction, valid=valid,
                           mask_reason=None if valid else 'insufficient_or_nonpositive_cell_or_control')
    return result


def closure(events, table):
    sums = defaultdict(lambda: {'gen': 0., 'raw': 0., 'single': 0., 'pair_supported': 0., 'masked_single': 0., 'masked_pairs': 0})
    vectors = []
    for event in events:
        if fold(event) != 1:
            continue
        weight = event['generator_weight'] * event['event_weight']
        local = defaultdict(lambda: {'gen': 0., 'raw': 0., 'single': 0., 'pair_supported': 0., 'masked_single': 0., 'masked_pairs': 0})
        for jet in event['jets']:
            for population in ('gen', 'reco'):
                for a, b in pairs(jet[population]):
                    base = base_key(event, jet, a, b)
                    dr = bin_index(separation(a, b), DR)
                    if base is None or dr is None:
                        continue
                    cell = table.get((*base, dr))
                    for n in (1, 2):
                        channel = (base[0], base[1], base[3], n, dr)  # activity integrated after weighted filling
                        row = local[channel]
                        if population == 'gen':
                            row['gen'] += weight * (a['pt'] * b['pt']) ** n
                        else:
                            row['raw'] += weight * (a['pt'] * b['pt']) ** n
                            single = pair_weight(a, b, n, weight)
                            row['single'] += single
                            if cell and cell['valid']:
                                row['pair_supported'] += single * cell['correction']
                            else:
                                row['masked_single'] += single
                                row['masked_pairs'] += 1
        for key, row in local.items():
            for name, value in row.items():
                sums[key][name] += value
        vectors.append({'id': list(physical_key(event)), 'fold': 1,
                        'channels': [{'key': list(k), **v} for k, v in sorted(local.items())]})
    spectra = []
    for key, row in sorted(sums.items()):
        complete = row['masked_pairs'] == 0
        spectra.append({'key': list(key), **row, 'pair': row['pair_supported'] if complete else None,
                        'pair_over_gen': row['pair_supported']/row['gen'] if complete and row['gen'] > 0 else None,
                        'status': 'NOT_VALIDATED' if complete else 'BLOCKED_INVALID_MAP_CELLS'})
    # Normalized-shape closure is separate from yields. Never normalize individual
    # background components; these are only truth/reco signal closure diagnostics.
    totals = defaultdict(lambda: {'gen': 0., 'raw': 0., 'single': 0., 'pair': 0., 'complete': True})
    for row in spectra:
        total = totals[tuple(row['key'][:-1])]
        for name in ('gen', 'raw', 'single'):
            total[name] += row[name]
        if row['pair'] is None:
            total['complete'] = False
        else:
            total['pair'] += row['pair']
    for row in spectra:
        total = totals[tuple(row['key'][:-1])]
        row['shape'] = {name: row[name]/total[name] if total[name] > 0 else None for name in ('gen', 'raw', 'single')}
        row['shape']['pair'] = row['pair']/total['pair'] if total['complete'] and total['pair'] > 0 else None
    return spectra, vectors


def run(events, min_pairs=10, control=CONTROL):
    if min_pairs < 1 or control[0] not in DR or control[1] not in DR or control[0] >= control[1]:
        raise ValueError('Control endpoints must be increasing stored dR edges; min_pairs positive')
    validate(events)
    table = derive(events, min_pairs, control)
    spectra, vectors = closure(events, table)
    return {'status': 'BLOCKED_NOT_A_PRODUCTION_PAYLOAD',
            'method': 'event-weighted n=0 single-track-corrected reco/common-reco-jet truth; C=alpha*gen/reco',
            'control': {'interval': list(control), 'status': 'PROPOSED_NOT_VALIDATED'},
            'min_raw_pairs_per_cell_and_control': min_pairs,
            'binning': {'pt': PT, 'jet': JET, 'activity_hiBin_provisional': ACTIVITY, 'dr_support': DR},
            'fold_algorithm': 'OOEventIdentity.h FNV1a64(namespace UTF8,zero,run/lumi/evt uint64 little endian) modulo 2; 0 train,1 holdout',
            'event_folds': [{'id': list(physical_key(e)), 'fold': fold(e)} for e in events],
            'map': [{'key': list(k), **v} for k, v in sorted(table.items())],
            'holdout_spectra': spectra, 'holdout_event_vectors': vectors,
            'covariance_status': 'NOT_RUN: common physical-event vectors retained for joint replicas; finite-MC map uncertainty not evaluated',
            'normalization_status': 'full stored-support shapes are diagnostic only; no final reporting interval chosen',
            'invalid_cell_policy': 'null correction and masked yields; no unity fallback, clamping, smoothing or negative-bin clipping',
            'physics_gates': ['single-track/reconstruction applicability', 'jet calibration applicability', 'large-angle control sensitivity', 'independent weighted closure all radii/powers', 'mixture/binning sensitivity', 'finite-MC covariance']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='JSONL events from export_forest.py')
    parser.add_argument('output', type=Path)
    parser.add_argument('--min-pairs', type=int, default=10)
    parser.add_argument('--control', nargs=2, type=float, default=CONTROL, metavar=('LOW', 'HIGH'))
    args = parser.parse_args()
    if args.min_pairs < 1:
        parser.error('--min-pairs must be positive')
    events = [json.loads(line) for line in args.input.read_text().splitlines() if line.strip()]
    result = run(events, args.min_pairs, tuple(args.control))
    result['tool_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['input_sha256'] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    result['input_path'] = str(args.input.resolve())
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'events': len(events), 'map_cells': len(result['map']), 'valid_map_cells': sum(c['valid'] for c in result['map']),
                      'holdout_channels': len(result['holdout_spectra']), 'status': result['status']}))


if __name__ == '__main__':
    main()
