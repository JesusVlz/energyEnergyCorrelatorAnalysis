import copy
import math
import unittest
from preflight import match, fold, physical_key, event_key


def events(n=12):
    return [{'sample': 'OO', 'sample_namespace': 'toy', 'id': [1, 1, i],
             'vz': 0., 'hf': 10., 'activity': 195, 'generator_weight': -1. if i % 3 == 0 else 1.,
             'jets': [{'radius': .4, 'pt': 50.}, {'radius': .8, 'pt': 90.}]} for i in range(n)]


class PoolTests(unittest.TestCase):
    def test_identity_fold_and_shared_radius_partners(self):
        rows = events()
        result = match(rows, 1., 20.)
        self.assertEqual(result['matched_signal_events'], len(rows))
        self.assertEqual(result['activity_coordinate']['status'], 'BLOCKED_DEGENERATE_COORDINATE')
        for row in result['assignments']:
            ids = [tuple(p['id']) for p in row['partners']]
            self.assertEqual(len(set(ids+[tuple(row['signal'])])), 3)
            self.assertTrue(all(p['fold'] == row['fold'] for p in row['partners']))
            self.assertEqual(len(row['jets']), 2)
        self.assertEqual(sum(r['uses'] for r in result['partner_reuse']), 2*len(rows))

    def test_independent_samples_do_not_mix(self):
        rows = events(1)
        for i in range(1, 12):
            other = copy.deepcopy(rows[0])
            other['sample_namespace'] = 'independent_'+str(i)
            rows.append(other)
        self.assertEqual(match(rows, 1., 20.)['matched_signal_events'], 0)

    def test_grid_matches_exhaustive_search_at_bin_boundaries(self):
        rows = events(40)
        for i, row in enumerate(rows):
            row['vz'] = (i % 7-3)*.49
            row['hf'] = (i % 9)*9.99
        result = match(rows, 1., 20.)
        found = {tuple(r['signal']): [tuple(p['id']) for p in r['partners']]
                 for r in result['assignments']}
        for signal in rows:
            candidates = [p for p in rows if physical_key(p) != physical_key(signal)
                          and fold(p) == fold(signal)
                          and abs(p['vz']-signal['vz']) <= 1.
                          and abs(p['hf']-signal['hf']) <= 20.]
            candidates.sort(key=lambda p: (abs(p['vz']-signal['vz'])+
                            abs(p['hf']-signal['hf'])/20., event_key(p), physical_key(p)))
            if len(candidates) >= 2:
                self.assertEqual(found[physical_key(signal)], [physical_key(p) for p in candidates[:2]])
            else:
                self.assertNotIn(physical_key(signal), found)

    def test_order_independent(self):
        self.assertEqual(match(events(), 1., 20.), match(list(reversed(events())), 1., 20.))

    def test_no_relaxation_for_sparse_pool(self):
        self.assertEqual(match(events(1), 1., 20.)['failed_signal_events'], 1)
        rows = events()
        for i, row in enumerate(rows):
            row['hf'] = i*100.
        self.assertEqual(match(rows, 1., 20.)['matched_signal_events'], 0)

    def test_duplicates_rejected(self):
        rows = events()
        with self.assertRaises(ValueError):
            match(rows+[copy.deepcopy(rows[0])], 1., 20.)

    def test_invalid_coordinates_and_tolerances(self):
        for bad in (0., -1., math.nan, math.inf):
            with self.assertRaises(ValueError):
                match(events(), bad, 20.)
        for key in ('vz', 'hf', 'generator_weight'):
            rows = events()
            rows[0][key] = math.nan
            with self.assertRaises(ValueError):
                match(rows, 1., 20.)

    def test_signed_weights_and_peripheral_not_clamped(self):
        result = match(events(), 1., 20.)
        self.assertTrue(any(r['signal_generator_weight'] < 0 for r in result['assignments']))
        self.assertTrue(all(r['signal_hiBin'] == 195 for r in result['assignments']))

    def test_zero_jet_events_remain_in_pool(self):
        rows = events()
        for row in rows[1:]:
            row['jets'] = []
        result = match(rows, 1., 20.)
        self.assertEqual(result['pool_events'], 12)
        self.assertEqual(result['matched_signal_events'], 1)


if __name__ == '__main__':
    unittest.main()
