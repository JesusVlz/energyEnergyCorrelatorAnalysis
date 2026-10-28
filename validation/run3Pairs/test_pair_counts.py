#!/usr/bin/env python3
"""Synthetic algorithm checks only; no production/physics-closure claim."""
import copy
import unittest
import export_forest as exporter
import pair_counts as p


def track(index, pt, phi, correction=1.):
    return {'index': index, 'pt': pt, 'eta': .2, 'phi': phi, 'correction': correction}


def jet(distance, reco=True, radius=.4):
    truth = [track(0, 3.1, 0), track(1, 4.1, distance)]
    return {'radius': radius, 'pt': 60., 'gen': copy.deepcopy(truth), 'reco': copy.deepcopy(truth) if reco else []}


def event(number, jets):
    return {'sample': 'OO', 'sample_namespace': 'OO-synthetic', 'id': [1, 2, number],
            'activity': 40, 'generator_weight': 3., 'event_weight': 2., 'jets': jets}


def event_in_fold(wanted, jets):
    return next(e for i in range(100) if p.fold(e := event(i, jets)) == wanted)


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
        self.assertEqual(small['gen'], 4 * 6.)  # generator/event factors each once
        self.assertEqual(small['single'], 2 * 6.)  # energy factor excluded
        self.assertEqual({row['key'][3] for row in result['holdout_spectra']}, {1, 2})
        for row in result['holdout_spectra']:
            self.assertAlmostEqual(row['pair_over_gen'], 1.)
        self.assertEqual(result['status'], 'BLOCKED_NOT_A_PRODUCTION_PAYLOAD')

    def test_holdout_never_enters_map(self):
        train = event_in_fold(0, [jet(.2)])
        holdout = event_in_fold(1, [jet(.01) for _ in range(20)])
        self.assertEqual(p.derive([train], 1), p.derive([train, holdout], 1))

    def test_missing_control_and_missing_cells_are_masked(self):
        train = event_in_fold(0, [jet(.01)])
        holdout = event_in_fold(1, [jet(.01), jet(.2)])
        result = p.run([train, holdout], 1)
        self.assertTrue(all(c['correction'] is None for c in result['map']))
        self.assertTrue(all(row['pair'] is None for row in result['holdout_spectra']))
        self.assertGreater(sum(row['masked_pairs'] for row in result['holdout_spectra']), 0)

    def test_radius_power_share_event_identity_and_holdout_vector(self):
        e = event_in_fold(1, [jet(.2), jet(.2, radius=.8)])
        before = p.event_key(e)
        e['jets'].reverse()
        self.assertEqual(before, p.event_key(e))
        _, vectors = p.closure([e], {})
        self.assertEqual(len(vectors), 1)
        self.assertEqual({tuple(row['key'][1:4:2]) for row in vectors[0]['channels']}, {('R04', 1), ('R04', 2), ('R08', 1), ('R08', 2)})
        other = copy.deepcopy(e)
        other['sample_namespace'] = 'independent-MC'
        self.assertNotEqual(p.event_key(e), p.event_key(other))

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

    def test_control_change_really_changes_candidate_mask(self):
        e = event_in_fold(0, [jet(.01), jet(.2)])
        self.assertTrue(all(row['valid'] for row in p.run([e], 1)['map']))
        self.assertTrue(all(not row['valid'] for row in p.run([e], 1, (.1, .15))['map']))

    def test_oo_lookup_cap_and_invalid_domain(self):
        class Axis:
            def GetXmin(self): return -3.
            def GetXmax(self): return 500.
        class Histogram:
            def __init__(self, value): self.value = value; self.coordinates = None
            def GetXaxis(self): return Axis()
            def GetYaxis(self): return Axis()
            def GetZaxis(self): return Axis()
            def GetDimension(self): return 3
            def FindBin(self, *values): self.coordinates = values; return 1
            def GetBinContent(self, unused): return self.value
        payload = [Histogram(v) for v in (.5, .1, .2, .1)]
        self.assertAlmostEqual(exporter.lookup(payload, 'OO', 50, .3, 198), 1.35)
        self.assertTrue(all(h.coordinates == (.3, 6.1, 198) for h in payload))
        with self.assertRaisesRegex(ValueError, 'hiBin'):
            exporter.lookup(payload, 'OO', 50, .3, 200)

    def test_r8_support_reaches_full_diameter_and_phi_wraps(self):
        self.assertIsNotNone(p.bin_index(1.5, p.DR))
        self.assertIsNone(p.bin_index(1.6, p.DR))
        self.assertAlmostEqual(p.separation(track(0, 3, math_pi_minus := 3.13), track(1, 4, -math_pi_minus)), 2*3.141592653589793-6.26)


if __name__ == '__main__':
    unittest.main()
