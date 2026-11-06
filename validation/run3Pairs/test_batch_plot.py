#!/usr/bin/env python3
"""Check batch moments against the established physical-event plotting path."""
import argparse
import copy
import json
import math
from pathlib import Path

import batch_plot
import plot_closure as plots


def fixture():
    dr = [0., .008, .02, .1, .4, .8]
    vectors = []
    for event_index in range(2):
        channels = []
        for n in (1, 2):
            truth = (1., 3. ** n)[event_index]
            raw = 2. ** n
            channels.append({'key': ['pp', 'R04', 0, n, 'ppRef', event_index + 1],
                             'gen': truth, 'raw': raw, 'single': raw * 1.1,
                             'pair_supported': raw * 1.2,
                             'pair_jet_independent_supported': raw * 1.15})
        vectors.append({'channels': channels})
    spectra = []
    for n in (1, 2):
        for angular in range(len(dr) - 1):
            key = ['pp', 'R04', 0, n, 'ppRef', angular]
            event_rows = [row for event in vectors for row in event['channels'] if row['key'] == key]
            row = {'key': key, 'status': 'SUPPORTED' if event_rows else 'EMPTY',
                   'masked_single': 0., 'masked_pairs': 0, 'reco_pairs': len(event_rows),
                   'fallback_pairs': 0}
            for field in plots.FIELDS:
                value = sum(x.get(plots.EVENT_FIELDS.get(field, field), 0.) for x in event_rows)
                row[field] = value if event_rows or field not in plots.EVENT_FIELDS else None
                if field in plots.EVENT_FIELDS:
                    row[plots.EVENT_FIELDS[field]] = value
            if n == 2 and angular == 1:
                row.update(pair=None, status='UNSUPPORTED', masked_single=.4, masked_pairs=1, reco_pairs=2)
            spectra.append(row)
    cell = {'valid': True, 'correction': .5, 'gen_pairs': 20, 'reco_pairs': 20,
            'large_angle_single_over_gen': 1.}
    return {'binning': {'dr_support': dr, 'jet': [40, 80], 'pt': [2, 4]},
            'map': [dict(cell, key=['pp', 'R04', 0, 0, 0, 0, 1])],
            'jet_independent_map': [dict(cell, key=['pp', 'R04', 0, 0, 0, 1])],
            'holdout_spectra': spectra, 'holdout_event_vectors': vectors,
            'jet_event_sums': [{'key': ['pp', 'R04', 0, 'ppRef'], 'weighted_jets': 2.}],
            'event_counts': {'raw': 5, 'selected': 5, 'train': 3, 'holdout': 2},
            'correction_version': 'run3-an751-v1', 'status': 'SYNTHETIC_EQUIVALENCE_ONLY'}


def moments_rows(statistics):
    return [dict(values, key=list(key)) for key, values in sorted(statistics.items())]


def numbers_equal(left, right):
    return (math.isnan(left) and math.isnan(right)) or math.isclose(left, right, rel_tol=1e-11, abs_tol=1e-12)


def compare_directory(left, right):
    compared = 0
    for key in left.GetListOfKeys():
        name = key.GetName()
        original = left.Get(name)
        candidate = right.Get(name)
        assert candidate, 'Missing ROOT object ' + name
        if original.InheritsFrom('TDirectory'):
            compared += compare_directory(original, candidate)
        elif original.InheritsFrom('TH1'):
            assert original.GetNbinsX() == candidate.GetNbinsX()
            for index in range(original.GetNbinsX() + 2):
                assert numbers_equal(original.GetBinContent(index), candidate.GetBinContent(index)), (name, index, 'content')
                assert numbers_equal(original.GetBinError(index), candidate.GetBinError(index)), (name, index, 'error')
            compared += 1
        elif original.InheritsFrom('TGraphErrors'):
            assert original.GetN() == candidate.GetN(), name
            for index in range(original.GetN()):
                for getter in ('GetPointX', 'GetPointY', 'GetErrorX', 'GetErrorY'):
                    assert numbers_equal(getattr(original, getter)(index), getattr(candidate, getter)(index)), (name, index, getter)
            compared += 1
    return compared


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    original = fixture()
    shape_result, shape_statistics, shape_support = plots.normalized_shapes(original)
    n1_first = ('pp', 'R04', 0, 1, 'ppRef', 1)
    error = plots.ratio_error(shape_statistics[n1_first], 'raw', .5, .25)
    assert math.isclose(error, math.sqrt(.5), rel_tol=1e-12), error
    assert shape_support[1]['suppressed_curves']['pair'] == 'unsupported populated bin in normalization interval'
    reduced = copy.deepcopy(original)
    reduced['event_statistics'] = moments_rows(plots.event_statistics(original))
    reduced['shape_statistics'] = moments_rows(shape_statistics)
    reduced['holdout_event_vectors'] = []
    reduced['event_folds'] = []
    absolute, normalized, supplied_shape, supplied_support = batch_plot.prepare_statistics(reduced)
    assert supplied_support == shape_support
    for key, values in shape_statistics.items():
        assert values == supplied_shape[key]
    bad = copy.deepcopy(reduced)
    bad['event_statistics'] = []
    try:
        batch_plot.prepare_statistics(bad)
    except ValueError as error:
        assert 'populated bin' in str(error)
    else:
        raise AssertionError('Missing populated-bin moments accepted')
    bad = copy.deepcopy(reduced)
    bad['event_statistics'][0]['gen2'] = -1.
    try:
        batch_plot.prepare_statistics(bad)
    except ValueError as error:
        assert 'negative second moment' in str(error)
    else:
        raise AssertionError('Negative second moment accepted')
    ROOT = plots.root_setup()
    reference_dir = args.output / 'reference'
    reference_dir.mkdir(exist_ok=True)
    reference = ROOT.TFile(str(reference_dir / 'closure.root'), 'RECREATE')
    plots.draw_closure(ROOT, original, reference_dir, reference)
    plots.draw_closure(ROOT, shape_result, reference_dir, reference, shape=True, statistics_override=shape_statistics)
    reference.Close()
    batch_dir = args.output / 'batch'
    report = batch_plot.render(ROOT, reduced, batch_dir, {'input': 'synthetic fixture'}, selected_jet=0)
    reference = ROOT.TFile.Open(str(reference_dir / 'closure.root'))
    candidate = ROOT.TFile.Open(str(batch_dir / 'closure.root'))
    compared = sum(compare_directory(reference.Get(name), candidate.Get(name)) for name in ('closure', 'closure_shape'))
    candidate.Close(); reference.Close()
    assert report['event_counts']['holdout'] == 2
    assert 'finite-training-map uncertainty remains pending' in report['uncertainties']
    summary = {'status': 'PASS', 'ROOT_version': ROOT.gROOT.GetVersion(),
               'compared_histograms_and_graphs': compared,
               'analytical_shape_ratio_error': math.sqrt(.5),
               'unsupported_pair_shape_suppressed': True,
               'missing_and_negative_moments_rejected': True,
               'core_plotter_sha256': batch_plot.file_sha256(Path(plots.__file__)),
               'adapter_sha256': batch_plot.file_sha256(Path(batch_plot.__file__))}
    (args.output / 'equivalence.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
