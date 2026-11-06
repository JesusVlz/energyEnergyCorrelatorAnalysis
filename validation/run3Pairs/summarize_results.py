#!/usr/bin/env python3
"""Summarize existing AN pair counts and final plots without adding diagnostics.

The pair integral is reported only when every populated angular bin in the
standard [.008,.4) range is supported. Supported contributions are explicitly
separate. Systems are loaded sequentially; event vectors are never retained
beyond the system being summarized.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gc
import hashlib
import json
import math
from pathlib import Path

DISPLAY = (.008, .4)
SELECTIONS = {'OO': ('MB', 'C0_30', 'C30_50', 'C50_90'), 'pp': ('ppRef',)}


def ratio(numerator, denominator):
    return numerator / denominator if numerator is not None and denominator > 0 else None


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        if rows:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def dump_json(value):
    return json.dumps(value, separators=(',', ':'), allow_nan=False)


def annotate_key(key, result, independent=False):
    """Physical bounds for an existing table key, without invented map cells."""
    pt, jets, activity, dr = (result['binning'][name] for name in ('pt', 'jet', 'activity_hiBin', 'dr_support'))
    has_dr = len(key) == (6 if independent else 7)
    normal = len(key) == (6 if independent else 7) or len(key) == (5 if independent else 6)
    a = key[2]
    j = None if independent else key[3]
    start = 3 if independent else 4
    high = key[start] if normal and isinstance(key[start], int) else None
    low = key[start+1] if high is not None and start+1 < len(key) else None
    ibin = key[-1] if has_dr else None
    return {'activity_bin': a,
            'activity_hiBin_low': activity[a] if key[0] == 'OO' and a is not None else None,
            'activity_hiBin_high': activity[a+1] if key[0] == 'OO' and a is not None else None,
            'jet_bin': j, 'jet_pt_low': jets[j] if j is not None else None,
            'jet_pt_high': jets[j+1] if j is not None else None,
            'pt_hi_bin': high, 'pt_hi_low': pt[high] if high is not None else None,
            'pt_hi_high': pt[high+1] if high is not None else None,
            'pt_lo_bin': low, 'pt_lo_low': pt[low] if low is not None else None,
            'pt_lo_high': pt[low+1] if low is not None else None,
            'dr_bin': ibin, 'dr_low': dr[ibin] if ibin is not None else None,
            'dr_high': dr[ibin+1] if ibin is not None else None}


def summarize_system(system, result, plots):
    dr, jets = result['binning']['dr_support'], result['binning']['jet']
    if DISPLAY[0] not in dr or DISPLAY[1] not in dr:
        raise ValueError('Stored angular edges do not support the exact AN display interval')
    display_bins = {i for i in range(len(dr)-1) if DISPLAY[0] <= dr[i] and dr[i+1] <= DISPLAY[1]}
    groups = defaultdict(list)
    for row in result['holdout_spectra']:
        if row['key'][0] == system and row['key'][-1] in display_bins:
            groups[tuple(row['key'][:-1])].append(row)
    jet_counts = {tuple(row['key']): row for row in result['jet_event_sums'] if row['key'][0] == system}
    plot_channels = {tuple(row['channel']): row for row in plots['closure']}
    shape_integrals = {tuple(row['channel']): row for row in plots['shape_support']}
    channels = []
    for selection in SELECTIONS[system]:
        for j in range(3):
            for n in (1, 2):
                key = (system, 'R04', j, n, selection)
                rows = groups[key]
                if len(rows) != len(display_bins) or key not in plot_channels or key not in shape_integrals:
                    raise ValueError('Missing final channel/bin/plot: '+str(key))
                totals = {field: sum(row[field] for row in rows) for field in
                          ('gen', 'raw', 'single', 'pair_supported', 'pair_jet_independent_supported',
                           'masked_single', 'masked_pairs', 'independent_masked_pairs',
                           'reco_pairs', 'gen_pairs', 'fallback_pairs', 'differential_pairs')}
                complete = totals['masked_pairs'] == 0
                independent_complete = totals['independent_masked_pairs'] == 0
                pair = totals['pair_supported'] if complete and totals['gen'] > 0 else None
                independent = totals['pair_jet_independent_supported'] if independent_complete and totals['gen'] > 0 else None
                # Cross-check the integrals already used by the standard plots.
                existing = shape_integrals[key]['integrals']
                for field, value in [('gen', totals['gen']), ('raw', totals['raw']), ('single', totals['single']), ('pair', pair)]:
                    normalized = value if value is not None and value > 0 else None
                    other = existing.get(field)
                    if normalized is None and other is None:
                        continue
                    if normalized is None or other is None or not math.isclose(normalized, other, rel_tol=1e-10, abs_tol=1e-10):
                        raise ValueError('Plot/count integral mismatch: '+str(key)+' '+field)
                exposure = jet_counts.get((system, 'R04', j, selection), {})
                channels.append({'system': system, 'selection': selection, 'radius': 'R04', 'jet_pt_low': jets[j],
                                 'jet_pt_high': jets[j+1], 'n': n, 'dr_low': DISPLAY[0], 'dr_high_exclusive': DISPLAY[1],
                                 'holdout_jets': exposure.get('jets', 0), 'weighted_holdout_jets': exposure.get('weighted_jets', 0.),
                                 'holdout_events_with_jets': exposure.get('events', 0),
                                 'truth_integral': totals['gen'], 'raw_over_truth': ratio(totals['raw'], totals['gen']),
                                 'single_over_truth': ratio(totals['single'], totals['gen']),
                                 'pair_over_truth': ratio(pair, totals['gen']),
                                 'pair_integral_complete': complete and totals['gen'] > 0,
                                 'pair_supported_contribution_over_truth': ratio(totals['pair_supported'], totals['gen']),
                                 'jet_independent_pair_over_truth': ratio(independent, totals['gen']),
                                 'supported_single_weight_fraction': ratio(totals['single']-totals['masked_single'], totals['single']),
                                 'fallback_pair_count': totals['fallback_pairs'], 'differential_pair_count': totals['differential_pairs'],
                                 'unsupported_pair_count': totals['masked_pairs'], 'reco_pair_count': totals['reco_pairs'],
                                 'gen_pair_count': totals['gen_pairs'],
                                 'unsupported_angular_bins': sum(row['masked_pairs'] > 0 for row in rows),
                                 'empty_angular_bins': sum(row['status'] == 'EMPTY' for row in rows),
                                 'zero_truth_angular_bins': sum(row['gen'] <= 0 for row in rows),
                                 'status': 'EMPTY' if totals['gen_pairs']+totals['reco_pairs'] == 0 else ('FULL_PAIR_SUPPORT' if complete else 'INCOMPLETE_PAIR_SUPPORT'),
                                 'closure_pdf': plot_channels[key]['pdf']})

    differential = {tuple(row['key']): row for row in result['map'] if row['key'][0] == system}
    independent = {tuple(row['key']): row for row in result['jet_independent_map'] if row['key'][0] == system}
    unsupported, consumed_controls = [], set()
    application_counts = Counter()
    for use in result['correction_application']:
        key = tuple(use['key'])
        if key[0] != system:
            continue
        source = use.get('correction_source') or 'unsupported'
        application_counts[source] += use['reco_pairs']
        cell = differential.get(key)
        indkey = (*key[:3], *key[4:]) if len(key) == 7 else None
        fallback = independent.get(indkey)
        if len(key) == 7:
            consumed_controls.add(('differential', *key[:-1]))
            consumed_controls.add(('jet_independent', *indkey[:-1]))
        if cell and cell['valid']:
            continue
        bounds = annotate_key(key, result) if len(key) == 7 else {name: None for name in annotate_key(('pp', 'R04', 0, 0, 0, 0, 0), result)}
        unsupported.append({'system': system, 'consumed_key': dump_json(key), **bounds,
                            'validation_reco_pair_count': use['reco_pairs'],
                            'validation_single_weight_n1': use['single_n1'], 'validation_single_weight_n2': use['single_n2'],
                            'effective_correction_source': source, 'effective_correction': use.get('correction'),
                            'differential_training_gen_pairs': cell['gen_pairs'] if cell else None,
                            'differential_training_reco_pairs': cell['reco_pairs'] if cell else None,
                            'differential_mask_reason': cell.get('mask_reason') if cell else 'MISSING_TRAINING_CELL',
                            'independent_training_gen_pairs': fallback['gen_pairs'] if fallback else None,
                            'independent_training_reco_pairs': fallback['reco_pairs'] if fallback else None,
                            'independent_mask_reason': fallback.get('mask_reason') if fallback else 'MISSING_TRAINING_CELL',
                            'fully_unsupported': source == 'unsupported'})

    controls, known_controls = [], set()
    angular_totals = defaultdict(lambda: Counter())
    for family, mapping in [('differential', differential), ('jet_independent', independent)]:
        for key, cell in mapping.items():
            total = angular_totals[(family, *key[:-1])]
            for field in ('gen_pairs', 'reco_pairs', 'unsupported_single_pairs'):
                total[field] += cell.get(field, 0)
    minimum = result['min_raw_pairs_per_cell_and_control']
    for control in result['control']['support']:
        key = tuple(control['key'])
        if key[0] != system:
            continue
        family = control['family']
        identity = (family, *key)
        known_controls.add(identity)
        reasons = []
        if control.get('unsupported_single_pairs', 0): reasons.append('control_unsupported_single_track_correction')
        for field in ('gen_pairs', 'reco_pairs'):
            if control[field] < minimum: reasons.append(f'control_{field}_below_{minimum}')
        for field in ('gen', 'single'):
            if control[field] <= 0: reasons.append(f'control_{field}_nonpositive')
        controls.append({'system': system, 'family': family, 'control_key': dump_json(key),
                         'interval_low_exclusive': .1, 'interval_high_exclusive': .4,
                         'control_gen_pair_count': control['gen_pairs'], 'control_reco_pair_count': control['reco_pairs'],
                         'control_gen_weight': control['gen'], 'control_raw_weight': control['raw'],
                         'control_single_weight': control['single'], 'control_single_over_gen': control.get('single_over_gen'),
                         'control_unsupported_single_pairs': control.get('unsupported_single_pairs', 0),
                         'angular_training_gen_pair_count': angular_totals[identity]['gen_pairs'],
                         'angular_training_reco_pair_count': angular_totals[identity]['reco_pairs'],
                         'consumed_by_validation_cells': identity in consumed_controls,
                         'control_supported': not reasons, 'mask_reason': ';'.join(reasons) if reasons else None})
    for identity in sorted(consumed_controls-known_controls):
        family, *key = identity
        controls.append({'system': system, 'family': family, 'control_key': dump_json(key),
                         'interval_low_exclusive': .1, 'interval_high_exclusive': .4,
                         'control_gen_pair_count': None, 'control_reco_pair_count': None,
                         'control_gen_weight': None, 'control_raw_weight': None, 'control_single_weight': None,
                         'control_single_over_gen': None, 'control_unsupported_single_pairs': None,
                         'angular_training_gen_pair_count': None, 'angular_training_reco_pair_count': None,
                         'consumed_by_validation_cells': True, 'control_supported': False,
                         'mask_reason': 'MISSING_TRAINING_CONTROL'})
    summary = {'system': system, 'correction_version': result['correction_version'], 'payload_status': result['status'],
               'channel_count': len(channels), 'fully_supported_channels_in_display_range': sum(row['pair_integral_complete'] for row in channels),
               'incomplete_channels_in_display_range': sum(row['status'] == 'INCOMPLETE_PAIR_SUPPORT' for row in channels),
               'input_selected_events': len(result['event_folds']),
               'training_events': result['derivation_coverage']['training_events'],
               'holdout_events': len(result['event_folds'])-result['derivation_coverage']['training_events'],
               'derivation_coverage': result['derivation_coverage'],
               'consumed_pair_sources_full_R4_support': dict(application_counts),
               'consumed_unsupported_differential_cells': len(unsupported),
               'consumed_cells_without_supported_fallback': sum(row['fully_unsupported'] for row in unsupported),
               'control_groups': len(controls), 'unsupported_control_groups': sum(not row['control_supported'] for row in controls),
               'map_support': {'differential_populated_cells': len(differential), 'differential_supported_cells': sum(row['valid'] for row in differential.values()),
                               'independent_populated_cells': len(independent), 'independent_supported_cells': sum(row['valid'] for row in independent.values())},
               'tables': plots['tables'],
               'closure_plots': len(plots['closure']), 'shape_plots': len(plots.get('shape_closure', [])),
               'resolution': plots.get('resolution', []),
               'observed_residuals_from_standard_plots': [{field: row[field] for field in ('channel', 'largest_single_track_residuals', 'largest_supported_pair_residuals')} for row in plots['closure']],
               'uncertainties': plots['uncertainties'], 'input_manifest': result.get('input_path')}
    return channels, unsupported, controls, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oo-pairs', type=Path, required=True)
    parser.add_argument('--pp-pairs', type=Path, required=True)
    parser.add_argument('--plots', type=Path, required=True, help='final directory containing OO/ and pp/ support_and_closure.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    inputs = [('OO', args.oo_pairs), ('pp', args.pp_pairs)]
    for system, path in inputs:
        if not path.is_file() or not (args.plots/system/'support_and_closure.json').is_file():
            parser.error('Final pair JSON and plot report required for '+system)
    channels, unsupported, controls, summaries = [], [], [], []
    for system, path in inputs:
        with path.open() as stream:
            result = json.load(stream)
        with (args.plots/system/'support_and_closure.json').open() as stream:
            plots = json.load(stream)
        system_channels, system_unsupported, system_controls, summary = summarize_system(system, result, plots)
        channels.extend(system_channels); unsupported.extend(system_unsupported); controls.extend(system_controls)
        summary['pair_json'] = str(path.resolve())
        summary['pair_json_sha256_from_plot_report'] = plots['input_sha256']
        summaries.append(summary)
        del result, plots, system_channels, system_unsupported, system_controls
        gc.collect()
    if len(channels) != 30:
        raise ValueError('Expected all 30 requested channels')
    args.output.mkdir(parents=True, exist_ok=True)
    for name, rows in [('channel_summary.csv', channels), ('unsupported_cells.csv', unsupported), ('control_support.csv', controls)]:
        write_csv(args.output/name, rows)
    summary = {'display_interval': list(DISPLAY), 'display_interval_endpoints': 'lower inclusive, upper exclusive',
               'normalization': 'ratio of integrated absolute momentum-weighted sums in common selected reconstructed jets',
               'pair_integral_policy': 'full pair ratio only if every populated display bin is supported; supported contribution is separately labeled and never a closure claim',
               'unsupported_cells_policy': 'actual validation-consumed invalid/missing differential cells; effective source distinguishes supported jet-independent fallback from unsupported',
               'control_policy': 'actual training counts including failures; missing training controls remain null',
               'channels': 30, 'systems': summaries,
               'outputs': [str((args.output/name).resolve()) for name in ('channel_summary.csv', 'unsupported_cells.csv', 'control_support.csv')],
               'tool_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output/'result_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'channels': len(channels), 'consumed_unsupported_differential_cells': len(unsupported),
                      'control_groups': len(controls), 'output': str(args.output.resolve())}))


if __name__ == '__main__':
    main()
