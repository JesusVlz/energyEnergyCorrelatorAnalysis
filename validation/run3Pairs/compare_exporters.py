#!/usr/bin/env python3
"""One bounded same-entry legacy/native reconciliation for AN pair inputs.

Use identical forest, namespace and entry range. Selection differences are
reported explicitly; numerical comparisons use common physical events/jets/
tracks only. This does not select or recalibrate a second analysis population.
"""
import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def read(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def identity(event):
    return (event['sample_namespace'], *event['id'])


def indexed(rows):
    result = {identity(row): row for row in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate physical event identity in bounded exporter comparison')
    return result


def summary(values):
    return {'count': len(values), 'min': min(values) if values else None,
            'max': max(values) if values else None, 'mean': math.fsum(values)/len(values) if values else None}


def compare(legacy_path, native_path):
    legacy_rows, native_rows = read(legacy_path), read(native_path)
    legacy, native = indexed(legacy_rows), indexed(native_rows)
    legacy_meta_path = Path(str(legacy_path)+'.metadata.json')
    native_meta_path = Path(str(native_path)+'.metadata.json')
    legacy_meta, native_meta = json.loads(legacy_meta_path.read_text()), json.loads(native_meta_path.read_text())
    if legacy_meta['input'] != native_meta['input'] or legacy_meta['sample_namespace'] != native_meta['sample_namespace']:
        raise ValueError('Bounded comparison requires identical forest path and namespace')
    if legacy_meta['entry_range'] != [0, native_meta['entry_limit_per_file']]:
        raise ValueError('Bounded comparison requires identical source-entry ranges')
    common = sorted(legacy.keys() & native.keys())
    counts = Counter(legacy_events=len(legacy), native_events=len(native), common_events=len(common),
                     legacy_only_events=len(legacy.keys()-native.keys()), native_only_events=len(native.keys()-legacy.keys()))
    values = {key: [] for key in ('legacy_generator_weight', 'native_generator_weight', 'native_stored_generator_weight',
                                  'legacy_total_weight', 'native_total_weight', 'native_over_legacy_total_weight',
                                  'native_over_stored_jet_pt', 'native_over_legacy_single_factor',
                                  'native_relative_eta_over_legacy_single_factor', 'source_event_hiBin_difference',
                                  'matched_track_pt_difference', 'matched_track_eta_difference', 'matched_track_phi_difference')}
    examples = {'event_weight_changes': [], 'jet_pt_changes': [], 'single_factor_changes': [], 'track_selection_changes': []}
    for key in common:
        old, new = legacy[key], native[key]
        if old['source_entry'] != new['source_entry']:
            raise ValueError('Same physical identity mapped to different source entries')
        for name in ('generator_weight',):
            values['legacy_'+name].append(old[name])
            values['native_'+name].append(new[name])
        values['native_stored_generator_weight'].append(new.get('stored_generator_weight', new['generator_weight']))
        values['source_event_hiBin_difference'].append(new['activity']-old['activity'])
        old_weight = old['generator_weight']*old['event_weight']
        new_weight = new['generator_weight']*new['event_weight']
        values['legacy_total_weight'].append(old_weight)
        values['native_total_weight'].append(new_weight)
        if old_weight:
            values['native_over_legacy_total_weight'].append(new_weight/old_weight)
        if not math.isclose(old_weight, new_weight, rel_tol=1e-6, abs_tol=1e-12):
            counts['changed_event_weights'] += 1
            if len(examples['event_weight_changes']) < 5:
                examples['event_weight_changes'].append({'id': list(key), 'source_entry': new['source_entry'],
                    'legacy_total': old_weight, 'native_total': new_weight,
                    'native_vz_weight': new.get('vz_weight'), 'native_centrality_weight': new.get('centrality_weight')})
        old_jets = {jet['index']: jet for jet in old['jets'] if math.isclose(jet['radius'], .4, rel_tol=0., abs_tol=1e-6)}
        new_jets = {jet['index']: jet for jet in new['jets'] if math.isclose(jet['radius'], .4, rel_tol=0., abs_tol=1e-6)}
        counts['legacy_jets_common_events'] += len(old_jets)
        counts['native_jets_common_events'] += len(new_jets)
        counts['common_jets'] += len(old_jets.keys() & new_jets.keys())
        counts['legacy_only_jets_common_events'] += len(old_jets.keys()-new_jets.keys())
        counts['native_only_jets_common_events'] += len(new_jets.keys()-old_jets.keys())
        for jet_index in old_jets.keys() & new_jets.keys():
            old_jet, new_jet = old_jets[jet_index], new_jets[jet_index]
            values['native_over_stored_jet_pt'].append(new_jet['pt']/old_jet['pt'])
            if not math.isclose(old_jet['pt'], new_jet['pt'], rel_tol=1e-6):
                counts['changed_jet_pt'] += 1
                if len(examples['jet_pt_changes']) < 5:
                    examples['jet_pt_changes'].append({'id': list(key), 'index': jet_index,
                        'legacy_stored_pt': old_jet['pt'], 'native_pt': new_jet['pt'], 'native_raw_pt': new_jet.get('raw_pt')})
            for population in ('reco', 'gen'):
                old_tracks = {track['index']: track for track in old_jet[population]}
                new_tracks = {track['index']: track for track in new_jet[population]}
                counts[population+'_legacy_tracks_common_jets'] += len(old_tracks)
                counts[population+'_native_tracks_common_jets'] += len(new_tracks)
                counts[population+'_common_tracks'] += len(old_tracks.keys() & new_tracks.keys())
                counts[population+'_legacy_only_tracks'] += len(old_tracks.keys()-new_tracks.keys())
                counts[population+'_native_only_tracks'] += len(new_tracks.keys()-old_tracks.keys())
                if old_tracks.keys() != new_tracks.keys() and len(examples['track_selection_changes']) < 5:
                    examples['track_selection_changes'].append({'id': list(key), 'jet_index': jet_index, 'population': population,
                        'legacy_only_indices': sorted(old_tracks.keys()-new_tracks.keys()), 'native_only_indices': sorted(new_tracks.keys()-old_tracks.keys())})
                for track_index in old_tracks.keys() & new_tracks.keys():
                    before, after = old_tracks[track_index], new_tracks[track_index]
                    for coordinate in ('pt', 'eta'):
                        values['matched_track_'+coordinate+'_difference'].append(after[coordinate]-before[coordinate])
                    values['matched_track_phi_difference'].append(math.remainder(after['phi']-before['phi'], 2*math.pi))
                    if population == 'reco':
                        ratio = after['correction']/before['correction']
                        values['native_over_legacy_single_factor'].append(ratio)
                        if 'legacy_relative_eta_correction' in after:
                            values['native_relative_eta_over_legacy_single_factor'].append(after['legacy_relative_eta_correction']/before['correction'])
                        if not math.isclose(ratio, 1., rel_tol=1e-6):
                            counts['changed_single_factors'] += 1
                            if len(examples['single_factor_changes']) < 5:
                                examples['single_factor_changes'].append({'id': list(key), 'jet_index': jet_index, 'track_index': track_index,
                                    'physical_eta': after['eta'], 'legacy_lookup_eta': before.get('single_track_lookup_eta'),
                                    'source_activity': new['activity'], 'legacy': before['correction'], 'native': after['correction']})
    ledger_path = Path(native_meta['native_outputs'])/'reco.jsonl.events.jsonl'
    ledger = {tuple(row['id']): row for row in read(ledger_path)} if ledger_path.exists() else {}
    removed = []
    for key in sorted(legacy.keys()-native.keys()):
        decision = ledger.get(key[1:], {})
        label = 'native_event_filters_fail' if decision.get('event_filters_pass') is False else 'other_native_selection_or_event_veto'
        counts['legacy_only_'+label] += 1
        removed.append({'id': list(key), 'source_entry': legacy[key]['source_entry'], 'reason': label})
    return {'status': 'BOUNDED_SAME_EVENT_RECONCILIATION', 'system': native_meta['system'],
            'input': native_meta['input'], 'sample_namespace': native_meta['sample_namespace'],
            'entry_range': legacy_meta['entry_range'], 'counts': dict(counts),
            'differences': {key: summary(value) for key, value in values.items()}, 'examples': examples,
            'legacy_only_events': removed,
            'native_only_events': [{'id': list(key), 'source_entry': native[key]['source_entry']} for key in sorted(native.keys()-legacy.keys())],
            'native_selection_weight_identity': native_meta['native_selection_weight_identity'],
            'native_eec_parity': native_meta['native_eec_parity'],
            'resolution': {'event_and_track_selection': 'analyzer-native selected events and common reco jets replace copied exporter cuts',
                           'event_weights': native_meta['event_weight_policy'], 'jet_calibration': native_meta['jet_policy'],
                           'single_track_correction': native_meta['single_track_coordinates']},
            'artifacts': {'legacy': str(legacy_path.resolve()), 'native': str(native_path.resolve()),
                          'legacy_metadata': str(legacy_meta_path.resolve()), 'native_metadata': str(native_meta_path.resolve())},
            'tool_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('legacy', type=Path)
    parser.add_argument('native', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = compare(args.legacy, args.native)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'system': result['system'], 'counts': result['counts'], 'native_eec_parity': result['native_eec_parity']}))


if __name__ == '__main__':
    main()
