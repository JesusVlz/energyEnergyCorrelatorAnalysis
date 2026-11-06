"""Check ledger exposure and failure gates without importing ROOT."""
import csv
from pathlib import Path
import tempfile
import unittest
from closure_statistics import JET_BINS, SELECTIONS
from plot_closure_pilot import read_ledgers


class ClosureLedgers(unittest.TestCase):
    def fixture(self, directory, status='matched', contribution=True, duplicate=False):
        prefix = Path(directory) / 'test'
        identity = 'OO:1:2:18446744073709551601'
        tables = {
            '.partners.csv': [dict(signal_key=identity, signal_background_key='18446744073709551599',
                                  signal_background_signature='4:0123456789abcdef:9876543210abcdef')],
            '.pool.csv': [dict(donor_key='OO:9:8:18446744073709551600', uses='1')],
            '.jets.csv': [dict(signal_key=identity, jet_index=0, jet_pt=60., centrality=10., event_weight=2., status=status)],
            '.contributions.csv': ([dict(signal_key=identity, jet_index=0, jet_pt=60., centrality=10., power=1,
                                   pairing=0, subevent=0, delta_r_bin=2, track_pt_bin=4, value=6.)] if contribution else [])}
        if duplicate:
            tables['.jets.csv'] *= 2
        for suffix, contents in tables.items():
            if not contents:
                fields = ['signal_key', 'jet_index', 'jet_pt', 'centrality', 'power', 'pairing',
                          'subevent', 'delta_r_bin', 'track_pt_bin', 'value']
            else:
                fields = list(contents[0])
            with Path(str(prefix)+suffix).open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader(); writer.writerows(contents)
        return prefix, identity

    def test_zero_cones_keep_matched_exposure(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix, _ = self.fixture(directory, contribution=False)
            result = read_ledgers(prefix, JET_BINS, SELECTIONS, 4)
            exposure = result['exposure']['selected_integrated', 0]
            self.assertEqual(exposure['selected_jets'], 1)
            self.assertEqual(exposure['matched_jets'], 1)
            self.assertEqual(exposure['matched_sumw'], 2.)
            self.assertFalse(result['event_vectors'])

    def test_failed_match_never_contributes_partial_estimator(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix, _ = self.fixture(directory, status='failed')
            with self.assertRaisesRegex(ValueError, 'outside common'):
                read_ledgers(prefix, JET_BINS, SELECTIONS, 4)

    def test_duplicate_jet_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix, _ = self.fixture(directory, duplicate=True)
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                read_ledgers(prefix, JET_BINS, SELECTIONS, 4)

    def test_large_ids_and_background_cluster_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix, identity = self.fixture(directory)
            result = read_ledgers(prefix, JET_BINS, SELECTIONS, 4)
            self.assertEqual(result['background_keys'][identity], '18446744073709551599')
            self.assertEqual(result['cluster_keys'][identity], 'background:18446744073709551599')
            vector = result['event_vectors']['selected_integrated', 0, 1, 2][identity]
            self.assertEqual(vector['PP'], 6.)
            # The output event weights are already in each analyzer contribution.
            # The reducer must not multiply by the jet's event weight again.
            self.assertEqual(vector['T'], 6.)

    def test_eec_threshold_independent_of_matching_multiplicity(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix, _ = self.fixture(directory)
            result = read_ledgers(prefix, JET_BINS, SELECTIONS, 5)
            self.assertFalse(result['event_vectors'])
            self.assertEqual(result['exposure']['selected_integrated', 0]['matched_jets'], 1)


if __name__ == '__main__':
    unittest.main()
