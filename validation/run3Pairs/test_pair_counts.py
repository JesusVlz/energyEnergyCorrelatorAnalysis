#!/usr/bin/env python3
"""Synthetic procedure checks, not a claim of Run-3 physics closure."""
import copy
import json
import unittest
import pair_counts as p


def track(index, pt, phi, correction=1.):
    return {'index': index, 'pt': pt, 'eta': .2, 'phi': phi, 'correction': correction}


def jet(distance, reco=True, radius=.4, pt=60.):
    truth = [track(0, 3.1, 0), track(1, 4.1, distance)]
    return {'radius': radius, 'pt': pt, 'gen': copy.deepcopy(truth), 'reco': copy.deepcopy(truth) if reco else []}


def event(number, jets):
    return {'sample': 'OO', 'sample_namespace': 'OO-synthetic', 'id': [1, 2, number],
            'activity': 40, 'generator_weight': 3., 'event_weight': 2., 'jets': jets}


def event_in_fold(wanted, jets):
    return next(e for i in range(100) if p.fold(e := event(i, jets)) == wanted)


def populated(result):
    return [row for row in result['holdout_spectra'] if row['status'] != 'EMPTY']


class PairCountsTest(unittest.TestCase):
    def test_shared_cpp_event_key_fixtures(self):
        for namespace, expected in [('OO2025_MC', 12337214294153387831), ('ppRef2024_MC', 15469188383063868161)]:
            self.assertEqual(p.event_key({'sample': 'OO', 'sample_namespace': namespace, 'id': [1, 2, 3]}), expected)

    def test_weights_once_and_independent_power(self):
        a, b = track(0, 3., 0, 2.), track(1, 4., .1, 4.)
        self.assertEqual(p.pair_weight(a, b, 1, 3., 5.), 1440.)
        self.assertEqual(p.pair_weight(a, b, 2, 3., 5.), 17280.)
        self.assertNotEqual(p.pair_weight(a, b, 2, 3., 5.), p.pair_weight(a, b, 1, 3., 5.)**2)

    def test_known_pair_loss_and_common_jet_holdout_closure(self):
        jets = [jet(.01, i < 2) for i in range(4)] + [jet(.2) for _ in range(4)]
        train, holdout = event_in_fold(0, jets), event_in_fold(1, jets)
        result = p.run([train, holdout], min_pairs=2)
        small = [c for c in result['map'] if c['key'][-1] == p.bin_index(.01, p.DR)][0]
        self.assertEqual(small['correction'], 2.)
        self.assertEqual(small['gen'], 4 * 6.)
        self.assertEqual(small['single'], 2 * 6.)
        self.assertEqual({row['key'][3] for row in populated(result)}, {1, 2})
        for row in populated(result):
            self.assertAlmostEqual(row['pair_over_gen'], 1.)
        self.assertEqual(result['status'], 'BLOCKED_NOT_A_PRODUCTION_PAYLOAD')
        json.dumps(result, allow_nan=False)

    def test_control_is_fixed_and_strict(self):
        train = event_in_fold(0, [jet(.01), jet(.1), jet(.4)])
        self.assertTrue(all(not c['valid'] for c in p.run([train], 1)['map']))
        train['jets'].append(jet(.2))
        result = p.run([train], 1)
        self.assertTrue(all(c['valid'] for c in result['map']))
        self.assertEqual(result['control']['support'][0]['gen_pairs'], 1)
        with self.assertRaisesRegex(ValueError, 'fixed'):
            p.run([train], 1, (.1, .15))

    def test_holdout_never_enters_either_map(self):
        train = event_in_fold(0, [jet(.2)])
        train.update(gen=jet(.2)['gen'], reco=jet(.2)['reco'])
        holdout = event_in_fold(1, [jet(.01) for _ in range(20)])
        holdout.update(gen=jet(.01)['gen'], reco=jet(.01)['reco'])
        self.assertEqual(p.derive_tables([train], 1), p.derive_tables([train, holdout], 1))

    def test_missing_control_and_missing_cells_are_masked(self):
        train = event_in_fold(0, [jet(.01)])
        holdout = event_in_fold(1, [jet(.01), jet(.2)])
        result = p.run([train, holdout], 1)
        self.assertTrue(all(c['correction'] is None for c in result['map']))
        self.assertTrue(all(row['pair'] is None for row in populated(result)))
        self.assertGreater(sum(row['masked_pairs'] for row in populated(result)), 0)
        self.assertTrue(all(row['correction'] is None for row in result['correction_application']))

    def test_supported_full_event_fallback_does_not_collapse_jets(self):
        train = event_in_fold(0, [])
        train.update(gen=jet(.2)['gen'], reco=jet(.2)['reco'])
        holdout = event_in_fold(1, [jet(.2, pt=100.)])
        result = p.run([train, holdout], 1)
        self.assertEqual(result['map'], [])
        self.assertEqual(len(result['jet_independent_map']), 1)
        for row in populated(result):
            self.assertEqual(row['fallback_pairs'], 1)
            self.assertEqual(row['pair_over_gen'], 1.)
            self.assertEqual(row['pair_jet_independent_over_gen'], 1.)
        self.assertEqual(result['correction_application'][0]['correction_source'], 'jet_independent_fallback')
        # Jet populations alone cannot fabricate the AN without-jets map.
        train.pop('gen')
        train.pop('reco')
        train['jets'] = [jet(.2)]
        result = p.run([train, holdout], 1)
        self.assertEqual(result['jet_independent_map'], [])
        self.assertTrue(all(row['pair'] is None for row in populated(result)))

    def test_an_independent_trigger_uses_jet_eta_acceptance(self):
        train = event_in_fold(0, [])
        low, high = track(0, 3.1, 0), track(1, 4.1, 0)
        low['eta'], high['eta'] = 1.5, 1.7
        train.update(gen=[low, high], reco=copy.deepcopy([low, high]))
        result = p.run([train], 1)
        self.assertFalse(result['jet_independent_map'])
        self.assertEqual(result['derivation_coverage']['independent_trigger_outside_jet_eta_pairs'], 2)
        for population in ('gen', 'reco'):
            train[population][0]['eta'], train[population][1]['eta'] = 1.7, 1.5
        result = p.run([train], 1)
        self.assertEqual(len(result['jet_independent_map']), 1)
        self.assertTrue(result['jet_independent_map'][0]['valid'])

    def test_unsupported_auxiliary_control_masks_entire_normalization_group(self):
        train = event_in_fold(0, [])
        gen = [track(0, 3.1, 0), track(1, 4.1, .01), track(2, 3.1, 1), track(3, 4.1, 1.2)]
        train.update(gen=gen, reco=copy.deepcopy(gen))
        train['reco'][3].update(correction=0., unsupported_single=True)
        result = p.run([train], 1)
        self.assertEqual(result['derivation_coverage']['input_unsupported_auxiliary_reco_tracks'], 1)
        self.assertEqual(result['derivation_coverage']['independent_unsupported_single_control_pairs'], 1)
        small = next(row for row in result['jet_independent_map'] if row['key'][-1] == p.bin_index(.01, p.DR))
        self.assertEqual(small['raw'], 6.)
        self.assertEqual(small['single'], 6.)
        self.assertTrue(small['single_is_complete'])
        self.assertFalse(small['valid'])
        self.assertIn('control_unsupported_single', small['mask_reason'])
        self.assertIsNone(small['large_angle_single_over_gen'])
        self.assertEqual(small['control']['raw'], 6.)
        self.assertEqual(small['control']['reco_pairs'], 1)
        self.assertEqual(small['control']['unsupported_single_raw'], 6.)
        self.assertFalse(small['control']['single_is_complete'])

    def test_unsupported_auxiliary_angular_cell_preserves_raw_without_poisoning_control(self):
        train = event_in_fold(0, [])
        gen = [track(0, 3.1, 0), track(1, 4.1, .01), track(2, 3.1, 1), track(3, 4.1, 1.2)]
        train.update(gen=gen, reco=copy.deepcopy(gen))
        train['reco'][1].update(correction=0., unsupported_single=True)
        result = p.run([train], 1)
        small = next(row for row in result['jet_independent_map'] if row['key'][-1] == p.bin_index(.01, p.DR))
        self.assertEqual(small['raw'], 6.)
        self.assertEqual(small['reco_pairs'], 1)
        self.assertEqual(small['single'], 0.)
        self.assertEqual(small['unsupported_single_pairs'], 1)
        self.assertFalse(small['single_is_complete'])
        self.assertFalse(small['valid'])
        self.assertEqual(small['large_angle_single_over_gen'], 1.)
        self.assertEqual(result['derivation_coverage']['independent_unsupported_single_control_pairs'], 0)
        self.assertTrue(any(row['valid'] for row in result['jet_independent_map']))

    def test_an_trigger_gate_can_exclude_unsupported_auxiliary_without_consuming_it(self):
        train = event_in_fold(0, [])
        low, high = track(0, 3.1, 0), track(1, 304., 0)
        low['eta'], high['eta'] = -1.5, -2.02
        train.update(gen=[low, high], reco=copy.deepcopy([low, high]))
        train['reco'][1].update(correction=0., unsupported_single=True)
        result = p.run([train], 1)
        self.assertFalse(result['jet_independent_map'])
        self.assertEqual(result['derivation_coverage']['input_unsupported_auxiliary_reco_tracks'], 1)
        self.assertEqual(result['derivation_coverage']['independent_unsupported_single_pairs'], 0)

    def test_unsupported_single_marker_never_allows_in_jet_or_unmarked_factor(self):
        train = event_in_fold(0, [jet(.2)])
        train['jets'][0]['reco'][0].update(correction=0., unsupported_single=True)
        with self.assertRaisesRegex(ValueError, 'only explicitly unsupported auxiliary'):
            p.run([train], 1)
        train['jets'] = []
        train.update(gen=jet(.2)['gen'], reco=jet(.2)['reco'])
        train['reco'][0]['correction'] = 0.
        with self.assertRaisesRegex(ValueError, 'Invalid single-track'):
            p.run([train], 1)

    def test_differential_preferred_and_factors_not_raised_to_n(self):
        train = event_in_fold(0, [jet(.01), jet(.01, False), jet(.2)])
        train.update(gen=jet(.01)['gen'] + [track(2, 3.1, 1), track(3, 4.1, 1.2)],
                     reco=jet(.01)['reco'] + [track(2, 3.1, 1), track(3, 4.1, 1.2)])
        holdout = event_in_fold(1, [jet(.01)])
        holdout['jets'][0]['reco'][0]['correction'] = 2.
        holdout['jets'][0]['reco'][1]['correction'] = 3.
        result = p.run([train, holdout], 1)
        for row in populated(result):
            self.assertEqual(row['differential_pairs'], 1)
            self.assertEqual(row['fallback_pairs'], 0)
            self.assertAlmostEqual(row['single_over_gen'], 6.)
            self.assertAlmostEqual(row['pair_over_gen'], 12.)

    def test_complete_scope_activity_boundaries_and_r4(self):
        events = []
        for activity in (0, 59, 60, 99, 100, 179, 180, 199):
            e = event_in_fold(1, [jet(.2), jet(.2, radius=.8)])
            e['id'][2] += 1000 * activity
            while p.fold(e) != 1:
                e['id'][2] += 1
            e['activity'] = activity
            events.append(e)
        result = p.run(events, 1)
        self.assertEqual(len(result['holdout_spectra']), 5 * 3 * 2 * (len(p.DR)-1))
        yields = {row['key'][4]: row['gen_pairs'] for row in populated(result) if row['key'][3] == 1}
        self.assertEqual(yields, {'MB': 8, 'C0_30': 2, 'C30_50': 2, 'C50_90': 2})
        self.assertTrue(all(row['key'][1] == 'R04' for row in result['holdout_spectra']))
        self.assertEqual(len(result['holdout_event_vectors']), 8)

    def test_native_float_radius_serialization(self):
        train = event_in_fold(0, [jet(.2)])
        holdout = event_in_fold(1, [jet(.2)])
        for source in (train, holdout):
            source['jets'][0]['radius'] = 0.4000000059604645
        result = p.run([train, holdout], 1)
        self.assertTrue(result['map'])
        self.assertTrue(all(row['pair_over_gen'] == 1. for row in populated(result)))

    def test_reject_duplicate_event_or_track(self):
        e = event(1, [jet(.2)])
        with self.assertRaisesRegex(ValueError, 'event identity'):
            p.validate([e, e])
        e['jets'][0]['reco'][1]['index'] = 0
        with self.assertRaisesRegex(ValueError, 'track index'):
            p.validate([e])

    def test_invalid_activity_not_silently_clamped(self):
        e = event(1, [jet(.2)])
        e['activity'] = 200
        with self.assertRaisesRegex(ValueError, 'Unsupported OO activity'):
            p.run([e], 1)

    def test_outside_table_pt_remains_in_truth_and_unsupported_reco(self):
        e = event_in_fold(1, [jet(.2)])
        for population in ('gen', 'reco'):
            e['jets'][0][population][0]['pt'] = 501.
        result = p.run([e], 1)
        for row in populated(result):
            self.assertGreater(row['gen'], 0)
            self.assertEqual(row['masked_pairs'], 1)
            self.assertIsNone(row['pair'])

    def test_nonpositive_weighted_cells_not_unity(self):
        train = event_in_fold(0, [jet(.2)])
        train['generator_weight'] = -1.
        result = p.run([train], 1)
        self.assertTrue(all(not row['valid'] and row['correction'] is None for row in result['map']))

    def test_physical_phi_wrap_and_r4_diameter(self):
        self.assertIsNotNone(p.bin_index(.799, p.DR))
        self.assertIsNone(p.bin_index(.8, p.DR))
        self.assertAlmostEqual(p.separation(track(0, 3, 3.13), track(1, 4, -3.13)), 2*3.141592653589793-6.26)


if __name__ == '__main__':
    unittest.main()
