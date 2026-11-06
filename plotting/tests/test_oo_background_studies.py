"""Numerical and selection regression checks that do not require PyROOT."""
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import Mock, patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plotOOMCBackgroundStudies as studies
import ooBackgroundDrawing as drawing


def example_histograms():
    """Two-bin toy sample with a 40-unit background and a 60-unit signal."""
    required = {key for terms in studies.COEFFICIENTS.values() for key in terms}
    raw = {key: (np.zeros(2), np.zeros(2)) for key in required}

    def store(key, value, variance=0.):
        raw[key] = (np.array([value, 2 * value], dtype=float),
                    np.array([variance, variance], dtype=float))

    store((0, 4), 100, 9)
    store((0, 0), 60)
    store((0, 1), 30)
    store((0, 3), 10)
    # The two cone methods use the same physical toy sample.
    for first, second, same1, same2, cross in ((3, 5, 4, 7, 6), (8, 10, 9, 12, 11)):
        store((first, 4), 40, 4)
        store((second, 4), 44, 16)
        store((same1, 4), 8, 36)
        store((same2, 4), 12, 64)
        store((cross, 4), 12, 121)
        store((first, 1), 21)
        store((second, 1), 23)
        store((first, 3), 19)
        store((second, 3), 21)
    return raw


class EstimatorTests(unittest.TestCase):
    def test_known_background_and_signal_for_both_methods(self):
        values, _ = studies.combine(example_histograms())
        np.testing.assert_allclose(values['Btrue'], [40, 80])
        for method in ('mix', 'perp'):
            np.testing.assert_allclose(values['B' + method], [40, 80])
            np.testing.assert_allclose(values['S' + method], [60, 120])

    def test_projected_signal_hijing_origin_is_already_merged(self):
        raw = example_histograms()
        # A separate S2 histogram is absent: projected S1 already has both orders.
        self.assertFalse(any(origin == 2 for _, origin in raw))
        values, _ = studies.combine(raw)
        np.testing.assert_allclose(values['PH'], [30, 60])
        np.testing.assert_allclose(values['PP'] + values['PH'] + values['HH'], values['T'])
        np.testing.assert_allclose(values['mix_PMavg'], [22, 44])

    def test_subtraction_keeps_negative_signal(self):
        raw = example_histograms()
        raw[0, 4] = (np.array([20., 40.]), np.array([9., 9.]))
        values, _ = studies.combine(raw)
        for name in ('Smix', 'Sperp'):
            np.testing.assert_allclose(values[name], [-20, -40])

    def test_variance_uses_squared_coefficients_including_subtractions(self):
        _, variances = studies.combine(example_histograms())
        # Four averaged terms contribute 30 and the subtracted cross term 121.
        for method in ('mix', 'perp'):
            np.testing.assert_allclose(variances['B' + method], [151, 151])
            np.testing.assert_allclose(variances['S' + method], [160, 160])


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.card = {
            'DataType': [10.], 'McCorrelationType': [3.],
            'BackgroundMethods': [6.], 'WeightExponent': [1.],
            'CentralityBinEdges': [-1., 0., 10., 90., 100.],
            'JetPtBinEdgesEEC': [60., 80.], 'TrackPtBinEdgesEEC': [1., 2.],
            'MinJetPtCut': [60.], 'MaxJetPtCut': [500.], 'JetRadius': [.4],
        }
        self.source = SimpleNamespace(IsZombie=lambda: False, Close=Mock())
        self.root = SimpleNamespace(TFile=SimpleNamespace(Open=Mock(return_value=self.source)))
        self.raw = example_histograms()
        self.reads = []

    def load(self, scope):
        def read_histogram(source, pairing, origin, cent, jet, track, njet):
            self.assertIs(source, self.source)
            self.assertNotEqual(origin, 2, 'Projected S1 must not be supplemented by S2')
            self.reads.append((pairing, origin, cent, jet, track))
            values, variances = self.raw[pairing, origin]
            # Distinct powers make omitted sentinel or peripheral slots detectable.
            factor = 10. ** cent
            return np.array([.008, .1, .39]), values * factor, variances * factor

        with patch.object(studies, 'read_card', return_value=self.card), \
                patch.object(studies, 'read_histogram', side_effect=read_histogram):
            return studies.load_channels(self.root, Path('toy.root'), scope)[1]

    def test_inclusive_folds_every_native_slot_without_a_centrality_selection(self):
        channels = self.load('integrated')
        self.assertEqual(len(channels), 2)  # one jet bin and the stored inclusive-jet bin
        for channel in channels:
            self.assertEqual(channel['centrality_indices'], [0, 1, 2, 3])
            self.assertEqual(channel['label'], 'Inclusive OO')
            np.testing.assert_allclose(channel['values']['T'], [111100, 222200])
            np.testing.assert_allclose(channel['values']['Smix'], [66660, 133320])
            np.testing.assert_allclose(channel['variances']['T'], [9999, 9999])
            self.assertNotIn('centrality', drawing._header(channel, 'closure').lower())
            self.assertEqual(drawing._group_label(channel, 'centrality'), 'Inclusive OO')
        self.assertEqual({entry[2] for entry in self.reads}, {0, 1, 2, 3})
        self.source.Close.assert_called_once()

    def test_split_keeps_native_slots_separate_including_sentinel_and_peripheral(self):
        channels = self.load('split')
        self.assertEqual(len(channels), 8)
        for cent in range(4):
            selected = [channel for channel in channels if channel['centrality_indices'] == [cent]]
            self.assertEqual(len(selected), 2)
            for channel in selected:
                np.testing.assert_allclose(channel['values']['T'], np.array([100, 200]) * 10. ** cent)
                self.assertIn('centrality', drawing._header(channel, 'closure').lower())
        self.assertEqual(channels[0]['centrality'], [-1., 0.])
        self.assertEqual(channels[-1]['centrality'], [90., 100.])
        self.source.Close.assert_called_once()


class RatioTests(unittest.TestCase):
    def test_invalid_denominators_and_nonfinite_inputs_are_masked(self):
        values, variances = drawing._ratio_arrays(
            [6., -4., 1., 1., np.inf, 1., np.nan],
            [3., 2., 0., -1., 2., np.nan, 2.],
            np.ones(7), np.ones(7))
        np.testing.assert_allclose(values[:2], [2., -2.])
        np.testing.assert_allclose(variances[:2], [5 / 9, 1.25])
        self.assertTrue(np.all(np.isnan(values[2:])))
        self.assertTrue(np.all(np.isnan(variances[2:])))


class FocusSelectionTests(unittest.TestCase):
    @staticmethod
    def channel(identifier, jet, threshold=2., inclusive=False):
        return dict(id=identifier, jet=jet, track_threshold=threshold,
                    jet_inclusive=inclusive)

    def setUp(self):
        self.selected = [self.channel('J{}-{}'.format(low, high), [low, high])
                         for low, high in ((40., 80.), (80., 120.),
                                           (120., 160.), (160., 200.))]
        self.excluded = [
            self.channel('no-J', [40., 200.], inclusive=True),
            self.channel('track-1', [40., 80.], threshold=1.),
            self.channel('track-3', [40., 80.], threshold=3.),
            self.channel('track-near-2', [40., 80.], threshold=2.000001),
            self.channel('below', [20., 40.]),
            self.channel('above', [200., 240.]),
            self.channel('overlap-below', [20., 60.]),
            self.channel('overlap-above', [160., 240.]),
        ]

    def test_default_focus_keeps_four_finite_jet_bins_at_track_threshold_two(self):
        channels = self.excluded[:3] + self.selected + self.excluded[3:]
        focused = studies.focus_channels(channels)
        self.assertEqual([c['id'] for c in focused], [c['id'] for c in self.selected])
        for actual, original in zip(focused, self.selected):
            self.assertIs(actual, original, 'Selection must preserve native channel metadata')

    def test_requested_range_and_threshold_are_used(self):
        alternate = self.channel('track-3', [80., 120.], threshold=3.)
        focused = studies.focus_channels(self.selected + [alternate], (80., 160.))
        self.assertEqual([c['id'] for c in focused], ['J80.0-120.0', 'J120.0-160.0'])
        self.assertEqual(studies.focus_channels(self.selected + [alternate],
                                               (80., 160.), 3.), [alternate])

    def test_jet_range_requires_increasing_boundaries(self):
        for jet_range in ((40., 40.), (200., 40.)):
            with self.subTest(jet_range=jet_range):
                with self.assertRaisesRegex(ValueError, 'increasing boundaries'):
                    studies.focus_channels(self.selected, jet_range)


class DrawingStyleTests(unittest.TestCase):
    def test_legend_has_exact_dataset_title_and_preserves_grid_free_axes(self):
        channel = dict(edges=np.array([.008, .02, .1, .4, .8]), radius=.4)
        figure, ax = drawing.plt.subplots()
        try:
            ax.plot([.02, .1], [1., 2.], label='True background')
            drawing._axes(ax, channel)
            drawing._legend(ax)
            figure.canvas.draw()
            legend = ax.get_legend()
            self.assertIsNotNone(legend)
            self.assertEqual(legend.get_title().get_text(), 'Pythia+Hijing OO 5.36 TeV')
            self.assertEqual([item.get_text() for item in legend.get_texts()], ['True background'])
            for axis in (ax.xaxis, ax.yaxis):
                ticks = axis.get_major_ticks() + axis.get_minor_ticks()
                self.assertFalse(any(tick.gridline.get_visible() for tick in ticks))
        finally:
            drawing.plt.close(figure)

    def test_axes_disable_major_and_minor_grids_on_both_axes(self):
        channel = dict(edges=np.array([.008, .02, .1, .4, .8]), radius=.4)
        with drawing.plt.rc_context({'axes.grid': True, 'axes.grid.which': 'both'}):
            figure, axes = drawing.plt.subplots(2, 1)
            try:
                for ax in axes:
                    ax.minorticks_on()
                    ax.grid(True, which='both', axis='both')
                figure.canvas.draw()
                self.assertTrue(any(line.get_visible() for ax in axes
                                    for line in ax.get_xgridlines() + ax.get_ygridlines()))
                for index, ax in enumerate(axes):
                    drawing._axes(ax, channel, xlabel=bool(index))
                figure.canvas.draw()
                for index, ax in enumerate(axes):
                    for axis in (ax.xaxis, ax.yaxis):
                        for tick_kind in ('major', 'minor'):
                            with self.subTest(panel=index, axis=axis.axis_name, ticks=tick_kind):
                                ticks = getattr(axis, 'get_{}_ticks'.format(tick_kind))()
                                self.assertTrue(ticks, 'Test requires actual rendered ticks')
                                self.assertFalse(any(tick.gridline.get_visible() for tick in ticks))
            finally:
                drawing.plt.close(figure)


if __name__ == '__main__':
    unittest.main()
