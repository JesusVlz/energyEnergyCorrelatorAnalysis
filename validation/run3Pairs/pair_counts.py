#!/usr/bin/env python3
"""AN Sec. 7.5.1 residual pair tables and independent-event R4 EEC closure.

The map uses n=0 pair counts. Detector factors multiply each EEC pair once,
independently of its momentum exponent. Unsupported corrections remain null.
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
ACTIVITY = (0., 60., 100., 180., 200.)  # source-event hiBin: 0-30/30-50/50-90/90-100%
DR = (0., .002, .005, .008, .01, .02, .04, .07, .1, .15, .2, .3, .4, .6, .8)
CONTROL = (.1, .4)
TRIGGER_ETA = 1.6  # native Run-3 R4 JetEtaCut; AN method 1 uses jet acceptance
VERSION = 'run3-an751-v2'
SELECTIONS = {'OO': ('MB', 'C0_30', 'C30_50', 'C50_90'), 'pp': ('ppRef',)}


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


def radius_is(jet, radius):
    return math.isclose(jet['radius'], radius, rel_tol=0., abs_tol=1e-6)


def activity_index(event):
    return bin_index(event['activity'], ACTIVITY) if event['sample'] == 'OO' else 0


def independent_key(event, a, b):
    pt_hi, pt_lo = sorted((a['pt'], b['pt']), reverse=True)
    indices = (activity_index(event), bin_index(pt_hi, PT), bin_index(pt_lo, PT))
    if any(i is None for i in indices):
        return None
    return (event['sample'], 'R04', *indices)


def base_key(event, jet, a, b):
    independent = independent_key(event, a, b)
    jet_bin = bin_index(jet['pt'], JET)
    if independent is None or jet_bin is None or not radius_is(jet, .4):
        return None
    return (*independent[:3], jet_bin, *independent[3:])


def selection_names(event):
    if event['sample'] == 'pp':
        return SELECTIONS['pp']
    activity = activity_index(event)
    return ('MB', SELECTIONS['OO'][activity + 1]) if activity < 3 else ('MB',)


def validate_population(tracks, allow_unsupported_auxiliary=False):
    for track in tracks:
        if not (math.isfinite(track['pt']) and 2 < track['pt'] and abs(track['eta']) < 2.4 and math.isfinite(track['phi'])):
            raise ValueError('Input violates strict measurement track acceptance')
        correction = track.get('correction', 1.)
        unsupported = track.get('unsupported_single', False)
        if unsupported:
            if not allow_unsupported_auxiliary or unsupported is not True or correction != 0:
                raise ValueError('Invalid single-track correction marker; only explicitly unsupported auxiliary reco may carry correction=0')
        elif not math.isfinite(correction) or correction <= 0:
            raise ValueError('Invalid single-track correction; no unity fallback')
    # Consume only index validation, not the O(N^2) combinations.
    indices = [track['index'] for track in tracks]
    if len(indices) != len(set(indices)):
        raise ValueError('Duplicate physical track index within one population')


def validate(events):
    seen = set()
    for event in events:
        identity = physical_key(event)
        if len(identity) != 4 or identity in seen:
            raise ValueError('Missing/duplicate physical event identity')
        seen.add(identity)
        if event['sample'] not in SELECTIONS or not identity[0] or any(v < 0 or v >= 1 << 64 for v in identity[1:]):
            raise ValueError('Unsupported system or invalid physical event identity')
        if event['sample'] == 'OO' and activity_index(event) is None:
            raise ValueError('Unsupported OO activity; no bin clamping or silent exclusion')
        for key in ('generator_weight', 'event_weight'):
            if not math.isfinite(event[key]):
                raise ValueError('Nonfinite event weight')
        if ('gen' in event) != ('reco' in event):
            raise ValueError('Jet-independent correction requires both event-level populations')
        if 'gen' in event:
            for population in ('gen', 'reco'):
                validate_population(event[population], allow_unsupported_auxiliary=(population == 'reco'))
        for jet in event['jets']:
            if bin_index(jet['pt'], JET) is None:
                raise ValueError('Unsupported reconstructed jet pT')
            if not any(radius_is(jet, value) for value in (.4, .8)):
                raise ValueError('Unsupported native radius')
            if not math.isfinite(jet.get('jet_weight', 1.)):
                raise ValueError('Nonfinite jet weight')
            for population in ('reco', 'gen'):
                validate_population(jet[population])


def count_cell():
    return {'gen': 0., 'single': 0., 'raw': 0., 'gen_sumw2': 0., 'single_sumw2': 0.,
            'gen_pairs': 0, 'reco_pairs': 0, 'unsupported_single_pairs': 0, 'unsupported_single_raw': 0.}


def fill_count(cell, population, a, b, weight):
    if population == 'gen':
        cell['gen'] += weight
        cell['gen_sumw2'] += weight * weight
        cell['gen_pairs'] += 1
    else:
        cell['raw'] += weight
        cell['reco_pairs'] += 1
        if a.get('unsupported_single', False) or b.get('unsupported_single', False):
            # Preserve the physical pair and its raw count. No corrected value
            # exists for this contribution: the containing cell is unsupported.
            cell['unsupported_single_pairs'] += 1
            cell['unsupported_single_raw'] += weight
        else:
            corrected = pair_weight(a, b, 0, weight)
            cell['single'] += corrected
            cell['single_sumw2'] += corrected * corrected


def finalize_table(cells, controls, min_pairs):
    result = {}
    for key, cell in sorted(cells.items()):
        control = controls.get(key[:-1], count_cell())
        alpha = control['single'] / control['gen'] if control['gen'] > 0 and not control['unsupported_single_pairs'] else None
        reasons = []
        for prefix, source in (('cell', cell), ('control', control)):
            if source['unsupported_single_pairs']:
                reasons.append(f'{prefix}_unsupported_single_track_correction')
            for name in ('gen_pairs', 'reco_pairs'):
                if source[name] < min_pairs:
                    reasons.append(f'{prefix}_{name}_below_{min_pairs}')
            for name in ('gen', 'single'):
                if source[name] <= 0:
                    reasons.append(f'{prefix}_{name}_nonpositive')
        valid = not reasons
        normalized_ratio = cell['single'] / cell['gen'] / alpha if valid else None
        result[key] = dict(cell, single_is_complete=not cell['unsupported_single_pairs'],
                           control=dict(control, single_is_complete=not control['unsupported_single_pairs']), large_angle_single_over_gen=alpha,
                           normalized_reco_over_gen=normalized_ratio,
                           correction=1. / normalized_ratio if valid else None, valid=valid,
                           mask_reason=None if valid else ';'.join(reasons))
    return result


def derive_tables(events, min_pairs=10, control=CONTROL):
    """Near-jet and full-event n=0 counts on the existing training-event fold.

    The jet-independent map uses selected source-event populations once per
    event with its higher-pT trigger inside jet eta acceptance (AN method 1). It is never synthesized by combining
    jet bins or deduplicating the union of jet constituents.
    """
    cells, controls = defaultdict(count_cell), defaultdict(count_cell)
    independent, independent_controls = defaultdict(count_cell), defaultdict(count_cell)
    coverage = {'training_events': 0, 'training_events_with_full_populations': 0,
                'differential_out_of_pt_domain_pairs': 0, 'independent_out_of_pt_domain_pairs': 0,
                'independent_trigger_outside_jet_eta_pairs': 0,
                'input_unsupported_auxiliary_reco_tracks': 0, 'training_unsupported_auxiliary_reco_tracks': 0,
                'independent_unsupported_single_pairs': 0, 'independent_unsupported_single_control_pairs': 0}
    for event in events:
        unsupported_tracks = sum(bool(track.get('unsupported_single', False)) for track in event.get('reco', []))
        coverage['input_unsupported_auxiliary_reco_tracks'] += unsupported_tracks
        if fold(event) != 0:
            continue
        coverage['training_unsupported_auxiliary_reco_tracks'] += unsupported_tracks
        coverage['training_events'] += 1
        weight = event['generator_weight'] * event['event_weight']
        populations = []
        for jet in event['jets']:
            if radius_is(jet, .4):
                populations.append((jet, cells, controls, weight * jet.get('jet_weight', 1.), 'differential'))
        if 'gen' in event and 'reco' in event:
            coverage['training_events_with_full_populations'] += 1
            populations.append((None, independent, independent_controls, weight, 'independent'))
        for jet, destination, normalization, local_weight, kind in populations:
            for population in ('gen', 'reco'):
                for a, b in pairs((event if jet is None else jet)[population]):
                    trigger = a if a['pt'] >= b['pt'] else b
                    if jet is None and abs(trigger['eta']) >= TRIGGER_ETA:
                        coverage['independent_trigger_outside_jet_eta_pairs'] += 1
                        continue
                    distance = separation(a, b)
                    dr = bin_index(distance, DR)
                    if dr is None:
                        continue
                    base = independent_key(event, a, b) if jet is None else base_key(event, jet, a, b)
                    if base is None:
                        coverage[f'{kind}_out_of_pt_domain_pairs'] += 1
                        continue
                    unsupported = population == 'reco' and (a.get('unsupported_single', False) or b.get('unsupported_single', False))
                    if unsupported:
                        coverage['independent_unsupported_single_pairs'] += 1
                    fill_count(destination[(*base, dr)], population, a, b, local_weight)
                    # Evaluate the actual pair coordinate: strict AN endpoints.
                    if control[0] < distance < control[1]:
                        if unsupported:
                            coverage['independent_unsupported_single_control_pairs'] += 1
                        fill_count(normalization[base], population, a, b, local_weight)
    return (finalize_table(cells, controls, min_pairs),
            finalize_table(independent, independent_controls, min_pairs), coverage)


def derive(events, min_pairs=10, control=CONTROL):
    """Compatibility API for the differential map; derive_tables returns both."""
    return derive_tables(events, min_pairs, control)[0]


def yield_cell():
    return {'gen': 0., 'raw': 0., 'single': 0., 'pair_supported': 0.,
            'pair_jet_independent_supported': 0., 'masked_single': 0., 'masked_pairs': 0,
            'independent_masked_single': 0., 'independent_masked_pairs': 0,
            'reco_pairs': 0, 'gen_pairs': 0, 'fallback_pairs': 0, 'differential_pairs': 0}


def correction_for(table, independent_table, base, dr):
    if base is None:
        return None, None, None
    cell = table.get((*base, dr))
    independent = independent_table.get((*base[:3], *base[4:], dr))
    independent_correction = independent['correction'] if independent and independent['valid'] else None
    if cell and cell['valid']:
        return cell['correction'], 'differential', independent_correction
    if independent_correction is not None:
        return independent_correction, 'jet_independent_fallback', independent_correction
    return None, None, None


def closure(events, table, independent_table=None, complete_grid=True, details=None):
    independent_table = independent_table or {}
    sums = defaultdict(yield_cell)
    jet_sums = defaultdict(lambda: {'weighted_jets': 0., 'jets': 0, 'events': 0})
    usage = defaultdict(lambda: {'reco_pairs': 0, 'single_n1': 0., 'single_n2': 0.})
    vectors = []
    if complete_grid:
        for sample, selections in SELECTIONS.items():
            for jet, n, selection, dr in itertools.product(range(len(JET)-1), (1, 2), selections, range(len(DR)-1)):
                sums[(sample, 'R04', jet, n, selection, dr)]
    for event in events:
        if fold(event) != 1:
            continue
        event_weight = event['generator_weight'] * event['event_weight']
        local, local_jets = defaultdict(yield_cell), defaultdict(lambda: {'weighted_jets': 0., 'jets': 0})
        for jet in event['jets']:
            if not radius_is(jet, .4):
                continue
            jet_bin = bin_index(jet['pt'], JET)
            weight = event_weight * jet.get('jet_weight', 1.)
            for selection in selection_names(event):
                counter = local_jets[(event['sample'], 'R04', jet_bin, selection)]
                counter['weighted_jets'] += weight
                counter['jets'] += 1
            for population in ('gen', 'reco'):
                for a, b in pairs(jet[population]):
                    base = base_key(event, jet, a, b)
                    dr = bin_index(separation(a, b), DR)
                    if dr is None:
                        raise ValueError('R4 pair lies outside the stored physical diameter')
                    correction, source, independent_correction = correction_for(table, independent_table, base, dr)
                    if population == 'reco':
                        usage_key = (*base, dr) if base is not None else (event['sample'], 'R04', activity_index(event), jet_bin, 'OUTSIDE_PT_DOMAIN', dr)
                        use = usage[usage_key]
                        use['reco_pairs'] += 1
                        use['single_n1'] += pair_weight(a, b, 1, weight)
                        use['single_n2'] += pair_weight(a, b, 2, weight)
                        use['correction_source'] = source
                        use['correction'] = correction
                    for n, selection in itertools.product((1, 2), selection_names(event)):
                        row = local[(event['sample'], 'R04', jet_bin, n, selection, dr)]
                        weighted = weight * (a['pt'] * b['pt']) ** n
                        if population == 'gen':
                            row['gen'] += weighted
                            row['gen_pairs'] += 1
                        else:
                            row['raw'] += weighted
                            row['reco_pairs'] += 1
                            single = pair_weight(a, b, n, weight)
                            row['single'] += single
                            if correction is not None:
                                row['pair_supported'] += single * correction
                                row['fallback_pairs' if source == 'jet_independent_fallback' else 'differential_pairs'] += 1
                            else:
                                row['masked_single'] += single
                                row['masked_pairs'] += 1
                            if independent_correction is not None:
                                row['pair_jet_independent_supported'] += single * independent_correction
                            else:
                                row['independent_masked_single'] += single
                                row['independent_masked_pairs'] += 1
        for key, row in local.items():
            for name, value in row.items():
                sums[key][name] += value
        for key, row in local_jets.items():
            for name, value in row.items():
                jet_sums[key][name] += value
            jet_sums[key]['events'] += 1
        vectors.append({'id': list(physical_key(event)), 'fold': 1,
                        'channels': [{'key': list(k), **v} for k, v in sorted(local.items())],
                        'jets': [{'key': list(k), **v} for k, v in sorted(local_jets.items())]})
    spectra = []
    for key, row in sorted(sums.items()):
        complete = row['masked_pairs'] == 0
        independent_complete = row['independent_masked_pairs'] == 0
        populated = row['reco_pairs'] > 0 or row['gen_pairs'] > 0
        result = {'key': list(key), **row,
                  'pair': row['pair_supported'] if complete and populated else None,
                  'pair_jet_independent': row['pair_jet_independent_supported'] if independent_complete and populated else None,
                  'supported_fraction_single': (row['single']-row['masked_single']) / row['single'] if row['single'] > 0 else None,
                  'status': 'EMPTY' if not populated else ('UNSUPPORTED_MAP_CELLS' if not complete else 'MEASURED_HOLDOUT')}
        for name in ('raw', 'single', 'pair', 'pair_jet_independent'):
            result[f'{name}_over_gen'] = result[name] / row['gen'] if result[name] is not None and row['gen'] > 0 else None
        spectra.append(result)
    totals = defaultdict(lambda: {'gen': 0., 'raw': 0., 'single': 0., 'pair': 0., 'complete': True})
    for row in spectra:
        total = totals[tuple(row['key'][:-1])]
        for name in ('gen', 'raw', 'single'):
            total[name] += row[name]
        if row['status'] != 'EMPTY':
            if row['pair'] is None:
                total['complete'] = False
            else:
                total['pair'] += row['pair']
    for row in spectra:
        total = totals[tuple(row['key'][:-1])]
        row['shape'] = {name: row[name]/total[name] if total[name] > 0 else None for name in ('gen', 'raw', 'single')}
        row['shape']['pair'] = row['pair']/total['pair'] if row['pair'] is not None and total['complete'] and total['pair'] > 0 else None
    if details is not None:
        details['correction_application'] = [{'key': list(k), **v} for k, v in sorted(usage.items(), key=lambda item: str(item[0]))]
    return spectra, vectors


def run(events, min_pairs=10, control=CONTROL):
    if min_pairs < 1 or tuple(control) != CONTROL:
        raise ValueError('AN normalization interval is fixed at 0.1 < DeltaR < 0.4; min_pairs positive')
    validate(events)
    table, independent, derivation_coverage = derive_tables(events, min_pairs, control)
    details = {}
    spectra, vectors = closure(events, table, independent, details=details)
    jet_sums = defaultdict(lambda: {'weighted_jets': 0., 'jets': 0, 'events': 0})
    for vector in vectors:
        for jet in vector['jets']:
            total = jet_sums[tuple(jet['key'])]
            total['events'] += 1
            for name in ('weighted_jets', 'jets'):
                total[name] += jet[name]
    control_rows = {}
    for family, mapping in (('differential', table), ('jet_independent', independent)):
        for key, cell in mapping.items():
            control_rows[(family, *key[:-1])] = {'family': family, 'key': list(key[:-1]),
                                                **cell['control'], 'single_over_gen': cell['large_angle_single_over_gen']}
    unsupported = [dict(key=row['key'], masked_pairs=row['masked_pairs'], masked_single=row['masked_single'],
                        single=row['single'], supported_fraction_single=row['supported_fraction_single'])
                   for row in spectra if row['masked_pairs']]
    return {'schema_version': 2, 'correction_version': VERSION,
            'status': 'BLOCKED_NOT_A_PRODUCTION_PAYLOAD',
            'method': 'AN751 n=0 single-track-corrected reco/common-reco-jet truth; C=alpha*gen/reco; alpha=sum(reco_single)/sum(gen) in strict control',
            'scope': {'radius': .4, 'systems': ['OO', 'ppRef'], 'track_pt_min_exclusive': 2., 'track_abs_eta_max_exclusive': 2.4,
                      'jet_pt': JET, 'powers': [1, 2], 'OO_selections': list(SELECTIONS['OO'])},
            'control': {'interval': list(control), 'endpoints': 'strict', 'status': 'FIXED_AN_0p1_0p4', 'support': list(control_rows.values())},
            'min_raw_pairs_per_cell_and_control': min_pairs,
            'binning': {'pt': PT, 'jet': JET, 'activity_hiBin': ACTIVITY, 'dr_support': DR},
            'map_key_fields': ['sample', 'radius', 'activity_bin', 'jet_bin', 'pt_hi_bin', 'pt_lo_bin', 'dr_bin'],
            'jet_independent_map_key_fields': ['sample', 'radius', 'activity_bin', 'pt_hi_bin', 'pt_lo_bin', 'dr_bin'],
            'holdout_key_fields': ['sample', 'radius', 'jet_bin', 'n', 'selection', 'dr_bin'],
            'fold_algorithm': 'OOEventIdentity.h FNV1a64(namespace UTF8,zero,run/lumi/evt uint64 little endian) modulo 2; 0 train,1 holdout',
            'event_folds': [{'id': list(physical_key(e)), 'fold': fold(e)} for e in events],
            'map': [{'key': list(k), **v} for k, v in sorted(table.items())],
            'jet_independent_map': [{'key': list(k), **v} for k, v in sorted(independent.items())],
            'jet_independent_population': 'AN method 1: full source-event selected gen/reco tracks; higher-pT trigger |eta|<1.6 (native JetEtaCut), associated |eta|<2.4; each pair once per event; no jet selection or jet-bin collapsing',
            'jet_independent_trigger_abs_eta_max_exclusive': TRIGGER_ETA,
            'AN_reference': '/eos/user/j/jvelazqu/HIN/Plan/AN2023_081_v13.pdf, Sec.7.5.1, methods1/2; existing TrackPairEfficiencyCorrector canonical trigger=high-pT member',
            'derivation_coverage': derivation_coverage,
            'holdout_spectra': spectra, 'holdout_event_vectors': vectors,
            'jet_event_sums': [{'key': list(k), **v} for k, v in sorted(jet_sums.items())],
            'unsupported_holdout_bins': unsupported,
            'correction_application': details['correction_application'],
            'covariance_status': 'common physical-event vectors retained; finite-training-MC map uncertainty not evaluated',
            'normalization_status': 'absolute momentum-weighted yields; optional full R4-support shapes reported separately; standard AN display uses 0.008 <= DeltaR < 0.4',
            'AN_shape_interval': [0.008, 0.4],
            'invalid_cell_policy': 'use supported full-event jet-independent cell when differential unsupported; otherwise null, never unity/clamp/smooth/clip',
            'unsupported_auxiliary_single_policy': 'explicit correction=0/unsupported_single=true allowed only in full-event reco auxiliary input; preserve raw counts, mask any eligible angular cell and normalization control consuming the pair; single sums in masked cells are supported contributions only; in-jet invalid factors fail',
            'background_contract': 'physical source-event track coordinates and corrections; residual map for same-source physical pairs only; never transplant to an artificial cross-event pair',
            'physics_gates': ['analyzer/exporter bounded same-event agreement', 'single-track/reconstruction applicability',
                              'jet calibration applicability', 'independent weighted closure in requested selections', 'training support']}


class JsonlEvents:
    """Re-iterable bounded-memory source for large native event exports."""
    def __init__(self, paths):
        self.paths = [paths] if isinstance(paths, Path) else list(paths)

    def __iter__(self):
        for path in self.paths:
            with path.open() as stream:
                for line in stream:
                    if line.strip():
                        yield json.loads(line)


def file_digest(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='JSONL events exported using analyzer selections and source-event corrections')
    parser.add_argument('output', type=Path)
    parser.add_argument('--input-list', action='store_true', help='input is a newline-separated JSONL path manifest; relative paths resolve beside the manifest')
    parser.add_argument('--min-pairs', type=int, default=10)
    parser.add_argument('--control', nargs=2, type=float, default=CONTROL, metavar=('LOW', 'HIGH'), help='fixed AN interval: 0.1 0.4')
    args = parser.parse_args()
    if args.min_pairs < 1 or tuple(args.control) != CONTROL:
        parser.error('--min-pairs must be positive; --control is fixed to 0.1 0.4')
    if args.input_list:
        paths = [Path(line.strip()) for line in args.input.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
        paths = [path if path.is_absolute() else args.input.parent/path for path in paths]
        if not paths:
            parser.error('input manifest is empty')
        if len({path.resolve() for path in paths}) != len(paths):
            parser.error('input manifest repeats a JSONL path')
    else:
        paths = [args.input]
    from preparation_contract import check_export_metadata
    tracking_contracts = check_export_metadata(paths)
    events = JsonlEvents(paths)
    result = run(events, args.min_pairs, tuple(args.control))
    result['tracking_contracts'] = tracking_contracts
    result['tool_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['input_sha256'] = file_digest(args.input)
    result['input_path'] = str(args.input.resolve())
    result['input_kind'] = 'JSONL_MANIFEST' if args.input_list else 'JSONL'
    result['input_files'] = [{'path': str(path.resolve()), 'sha256': file_digest(path)} for path in paths]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'events': len(result['event_folds']), 'map_cells': len(result['map']), 'valid_map_cells': sum(c['valid'] for c in result['map']),
                      'jet_independent_cells': len(result['jet_independent_map']),
                      'valid_jet_independent_cells': sum(c['valid'] for c in result['jet_independent_map']),
                      'holdout_channels': len(result['holdout_spectra']), 'unsupported_holdout_bins': len(result['unsupported_holdout_bins']),
                      'status': result['status']}))


if __name__ == '__main__':
    main()
