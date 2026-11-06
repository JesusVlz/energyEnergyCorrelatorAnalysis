"""Synthetic manager regression: no mixed-cone shape stabilization in raw view."""
from array import array
from pathlib import Path
import tempfile
import unittest

from closure_statistics import ALL, JET_BINS, SELECTIONS
from plot_closure_pilot import initialize_root, project


class RawManagerProjection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ROOT = initialize_root()

    def test_two_powers_origin_merge_and_raw_jet_bin_shapes(self):
        ROOT = self.ROOT
        with tempfile.TemporaryDirectory(prefix='eec-manager-test-') as directory:
            path = str(Path(directory) / 'synthetic.root')
            file = ROOT.TFile.Open(path, 'CREATE')
            card = file.mkdir('JCard'); card.cd()
            metadata = dict(DataType=[10.], McCorrelationType=[3.], JetRadius=[.4],
                BackgroundMethods=[6.], WeightExponent=[1., 2.],
                CentralityBinEdges=[-1., 0., 30., 50., 90., 100.],
                TrackPtBinEdges=[.7, 300.], JetPtBinEdgesEEC=[40., 80., 120., 160.],
                TrackPtBinEdgesEEC=[2., 3., 300.])
            for name, values in metadata.items():
                vector = ROOT.TVectorF(1, len(values), array('f', values))
                vector.Write(name)
            file.cd()
            edges = ([0., .008, .05, .39, .8], [40., 80., 120., 160.], [2., 3., 300.],
                     [-1., 0., 30., 50., 90., 100.],
                     [i - .5 for i in range(14)], [i - .5 for i in range(5)], [-.5, .5, 1.5])
            sparse = ROOT.THnSparseD('energyEnergyCorrelator', '', 7,
                array('i', [len(value)-1 for value in edges]),
                array('d', [value[0] for value in edges]),
                array('d', [value[-1] for value in edges]))
            sparse.Sumw2()
            for axis, values in enumerate(edges):
                sparse.GetAxis(axis).Set(len(values)-1, array('d', values))
            for power in (1, 2):
                for origin, value in ((0, 2.), (1, 3.), (2, 5.), (3, 7.)):
                    sparse.Fill(array('d', [.02, 60., 2.5, 10., 0., origin, power - 1.]), value * power)
                # Different per-jet-bin MM shapes: the ordinary manager load
                # would mix these shapes. The raw closure entry point must not.
                sparse.Fill(array('d', [.02, 60., 2.5, 10., 4., 3., power - 1.]), 2. * power)
                sparse.Fill(array('d', [.2, 100., 2.5, 10., 4., 3., power - 1.]), 8. * power)
            sparse.Write(); file.Close()
            source = ROOT.TFile.Open(path, 'READ')
            projected, retained = project(ROOT, source, JET_BINS, SELECTIONS, (1, 2), 2.)
            self.assertEqual(len(projected), 24)
            for power in (1, 2):
                categories, metrics = projected['C0_30', 0, power]
                first = metrics['T'].FindFixBin(.02)
                second = metrics['T'].FindFixBin(.2)
                self.assertAlmostEqual(metrics['PH'].GetBinContent(first) * metrics['PH'].GetBinWidth(first), 8. * power)
                self.assertAlmostEqual(metrics['T'].GetBinContent(first) * metrics['T'].GetBinWidth(first), 17. * power)
                self.assertAlmostEqual(metrics['truth_sum_residual'].GetBinContent(first), 0.)
                self.assertAlmostEqual(categories[4, ALL].GetBinContent(first) * categories[4, ALL].GetBinWidth(first), 2. * power)
                self.assertEqual(categories[4, ALL].GetBinContent(second), 0.)
                self.assertEqual(projected['C30_50', 0, power][1]['T'].Integral(), 0.)
                self.assertAlmostEqual(projected['selected_integrated', 0, power][1]['T'].Integral('width'), 17. * power)
            self.assertEqual(ROOT.EECCard(source).GetNWeightExponents(), 2)
            # Keep managers/cards alive until all projected histograms are read.
            self.assertEqual(len(retained), 4)


if __name__ == '__main__':
    unittest.main()
