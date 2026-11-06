#!/usr/bin/env python3
"""Reproducible diagnostics of merged OO jet-pT / EEC response projections.

ROOT response X is reco and Y is truth. Arrays retain that order, i.e.
response[reco, truth]; figures explicitly display truth horizontally.
"""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil

import numpy as np


def safe_divide(numerator, denominator):
    numerator, denominator = np.broadcast_arrays(
        np.asarray(numerator, dtype=float), np.asarray(denominator, dtype=float))
    return np.divide(numerator, denominator, out=np.full(numerator.shape, np.nan),
                     where=np.isfinite(denominator) & (denominator != 0))


def migration_metrics(matrix, reported_bins):
    """Matched-only metrics, including all stored guard bins in denominators."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError('Response must be square')
    if not np.isfinite(matrix).all() or np.any(matrix < 0):
        raise ValueError('These probability diagnostics require finite nonnegative weights')
    reported = np.asarray(reported_bins, dtype=int)
    if not np.array_equal(reported, np.arange(reported[0], reported[-1] + 1)):
        raise ValueError('Reported bins must be contiguous')
    rows, cols = matrix.sum(axis=1), matrix.sum(axis=0)
    groups = [np.arange(reported[0]), reported,
              np.arange(reported[-1] + 1, matrix.shape[0])]
    feed_in = np.array([[matrix[r, g].sum() for g in groups] for r in reported])
    feed_out = np.array([[matrix[g, t].sum() for g in groups] for t in reported])
    return {
        'conditional': safe_divide(matrix, cols[None, :]),
        'purity': safe_divide(np.diag(matrix), rows),
        'stability': safe_divide(np.diag(matrix), cols),
        'feed_in': safe_divide(feed_in, rows[reported, None]),
        'feed_out': safe_divide(feed_out, cols[reported, None]),
    }


def collapse_blocks(matrix, n_jet, n_dr):
    matrix = np.asarray(matrix)
    if matrix.shape != (n_jet * n_dr, n_jet * n_dr):
        raise ValueError('Flattened response dimensions do not match jet/angular binning')
    return matrix.reshape(n_jet, n_dr, n_jet, n_dr).sum(axis=(1, 3))


def moment_correlation(covariance):
    """Diagonal normalization of the stored within-jet second moments."""
    covariance = np.asarray(covariance, dtype=float)
    if not np.isfinite(covariance).all() or np.any(np.diag(covariance) < 0):
        raise ValueError('Invalid covariance entries')
    if not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-10):
        raise ValueError('Stored covariance is not symmetric')
    denominator = np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    correlation = safe_divide(covariance, denominator)
    if np.any(np.abs(correlation[np.isfinite(correlation)]) > 1 + 1e-7):
        raise ValueError('Stored second moments violate correlation bounds')
    return correlation


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def axis_edges(axis):
    return np.array([axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)])


def card_vector(root_file, name):
    obj = root_file.Get('JCard/' + name)
    if not obj:
        raise ValueError('Missing JCard/' + name)
    return np.array([obj.GetMatrixArray()[i] for i in range(obj.GetNoElements())])


def read_histogram(root_file, path, audits):
    hist = root_file.Get(path)
    if not hist or not hist.InheritsFrom('TH1'):
        raise ValueError('Missing histogram: ' + path)
    nx = hist.GetNbinsX()
    if hist.GetDimension() == 2:
        ny = hist.GetNbinsY()
        all_bins = np.array([[hist.GetBinContent(i, j) for j in range(ny + 2)]
                             for i in range(nx + 2)])
        values = all_bins[1:-1, 1:-1].copy()
        errors = np.array([[hist.GetBinError(i, j) for j in range(1, ny + 1)]
                           for i in range(1, nx + 1)])
        flow = np.concatenate([all_bins[0, :], all_bins[-1, :],
                               all_bins[1:-1, 0], all_bins[1:-1, -1]])
        y_edges = axis_edges(hist.GetYaxis())
    else:
        all_bins = np.array([hist.GetBinContent(i) for i in range(nx + 2)])
        values = all_bins[1:-1].copy()
        errors = np.array([hist.GetBinError(i) for i in range(1, nx + 1)])
        flow, y_edges = all_bins[[0, -1]], None
    if not np.isfinite(all_bins).all() or np.any(values < 0):
        raise ValueError('Nonfinite or negative histogram contents in ' + path)
    audit = {'path': path, 'shape': list(values.shape), 'regular_sum': float(values.sum()),
             'flow_sum': float(flow.sum()), 'absolute_flow_sum': float(abs(flow).sum()),
             'x_edges': axis_edges(hist.GetXaxis()).tolist()}
    if y_edges is not None:
        audit['y_edges'] = y_edges.tolist()
    audits.append(audit)
    if np.any(flow != 0):
        raise ValueError('Nonzero under/overflow requires an explicit plotting policy: ' + path)
    return values, errors, axis_edges(hist.GetXaxis()), y_edges


def load_inputs(input_dir, track_pt):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    receipt_path = input_dir / 'finalization_receipt.json'
    receipt = json.loads(receipt_path.read_text())
    data = {'folds': ['train', 'validation'], 'powers': [1, 2],
            'jet': {}, 'eec': {}, 'metrics': {}}
    provenance = {'files': [], 'receipt_sha256': sha256(receipt_path),
                  'manifest_sha256': receipt['manifest_sha256'],
                  'merged_archives': receipt['files'], 'raw_entries': receipt['raw_entries'],
                  'root_version': ROOT.gROOT.GetVersion()}
    for fold in data['folds']:
        data['eec'][fold], data['metrics'][fold] = {}, {}
        for power in data['powers']:
            path = input_dir / ('OO_%s_inclusive_N%d.root' % (fold, power))
            digest = sha256(path)
            if digest != receipt['outputs'][path.name]['sha256']:
                raise ValueError('Input checksum differs from finalization receipt: ' + str(path))
            source = {'path': str(path.absolute()), 'sha256': digest, 'histograms': []}
            provenance['files'].append(source)
            root_file = ROOT.TFile.Open(str(path))
            if not root_file or root_file.IsZombie():
                raise ValueError('Cannot open ' + str(path))
            try:
                jets = card_vector(root_file, 'JetPtBinEdgesEEC')
                thresholds = card_vector(root_file, 'TrackPtBinEdgesEEC')
                track_indices = np.flatnonzero(np.isclose(thresholds[:-1], track_pt))
                if len(track_indices) != 1:
                    raise ValueError('Track threshold must match an embedded card bin edge')
                track_index = int(track_indices[0])
                first = int(card_vector(root_file, 'FirstUnfoldedJetPtBin')[0])
                last = int(card_vector(root_file, 'LastUnfoldedJetPtBin')[0])
                source['card'] = {name: card_vector(root_file, name).tolist() for name in [
                    'JetRadius', 'WeightExponent', 'MinJetPtUnfoldingReco',
                    'MinJetPtUnfoldingTruth', 'CentralityBinEdges', 'MinTrackPtCut',
                    'MinPairTrackPt', 'JetPtBinEdgesEEC', 'TrackPtBinEdgesEEC',
                    'FirstUnfoldedJetPtBin', 'LastUnfoldedJetPtBin']}
                source['track_threshold_index'] = track_index
                if len(card_vector(root_file, 'CentralityBinEdges')) != 2:
                    raise ValueError('Expected inclusive projection')
                if 'jet_edges' not in data:
                    data['jet_edges'] = jets
                    data['reported_bins'] = list(range(first, last + 1))
                np.testing.assert_array_equal(data['jet_edges'], jets)
                if data['reported_bins'] != list(range(first, last + 1)):
                    raise ValueError('Inconsistent reported jet bins')
                jet, eec = {}, {}
                for stem, destination, suffix in [
                    ('oneDimensionalJetPtUnfolding', jet, '_C0'),
                    ('jetPtUnfolding', eec, '_C0T%d' % track_index)]:
                    kinds = ['Measured', 'Truth', 'Response']
                    if destination is eec:
                        kinds.append('Covariance')
                    for kind in kinds:
                        name = stem + kind
                        values, errors, x_edges, y_edges = read_histogram(
                            root_file, name + '/' + name + suffix, source['histograms'])
                        destination[kind.lower()] = values
                        destination[kind.lower() + '_sumw2'] = errors ** 2
                        if destination is jet:
                            np.testing.assert_allclose(x_edges, jets, rtol=0, atol=1e-10)
                            if y_edges is not None:
                                np.testing.assert_allclose(y_edges, jets, rtol=0, atol=1e-10)
                        else:
                            if 'flattened_edges' not in data:
                                n_jet = len(jets) - 1
                                if (len(x_edges) - 1) % n_jet:
                                    raise ValueError('Nonintegral angular blocks')
                                n_dr = (len(x_edges) - 1) // n_jet
                                data['flattened_edges'] = x_edges
                                data['dr_edges'] = x_edges[:n_dr + 1].copy()
                                for block in range(n_jet):
                                    np.testing.assert_allclose(
                                        x_edges[block * n_dr:block * n_dr + n_dr + 1],
                                        data['dr_edges'] + block * data['dr_edges'][-1],
                                        rtol=0, atol=1e-9)
                            np.testing.assert_allclose(x_edges, data['flattened_edges'], rtol=0, atol=1e-10)
                            if y_edges is not None:
                                np.testing.assert_allclose(y_edges, data['flattened_edges'], rtol=0, atol=1e-10)
                if fold in data['jet']:
                    for key in jet:
                        np.testing.assert_array_equal(data['jet'][fold][key], jet[key])
                else:
                    data['jet'][fold] = jet
                data['eec'][fold][power] = eec
            finally:
                root_file.Close()
        data['metrics'][fold]['jet'] = migration_metrics(data['jet'][fold]['response'], data['reported_bins'])
        for power in data['powers']:
            eec = data['eec'][fold][power]
            collapsed = collapse_blocks(eec['response'], len(data['jet_edges']) - 1,
                                        len(data['dr_edges']) - 1)
            metrics = migration_metrics(collapsed, data['reported_bins'])
            metrics['jet_migration'] = collapsed
            metrics['conditional'] = safe_divide(eec['response'], eec['response'].sum(axis=0)[None, :])
            metrics['correlation'] = moment_correlation(eec['covariance'])
            data['metrics'][fold][power] = metrics
    return data, provenance


def validate_metrics(data):
    checks = []
    for fold in data['folds']:
        for channel in ['jet'] + data['powers']:
            metric = data['metrics'][fold][channel]
            conditional = metric['conditional']
            supported = np.isfinite(conditional).any(axis=0)
            np.testing.assert_allclose(conditional[:, supported].sum(axis=0), 1., rtol=1e-12)
            for key in ['feed_in', 'feed_out']:
                array = metric[key]
                support = np.isfinite(array).any(axis=1)
                np.testing.assert_allclose(array[support].sum(axis=1), 1., rtol=1e-12)
            for key in ['conditional', 'purity', 'stability', 'feed_in', 'feed_out']:
                finite = metric[key][np.isfinite(metric[key])]
                if np.any(finite < 0) or np.any(finite > 1 + 1e-12):
                    raise ValueError('Metric outside probability bounds')
            checks.append({'fold': fold, 'channel': channel,
                           'supported_conditional_columns': int(supported.sum()),
                           'total_conditional_columns': len(supported), 'status': 'PASS'})
    return checks


def flatten_arrays(value, prefix='', result=None):
    if result is None:
        result = {}
    if isinstance(value, dict):
        for key, item in value.items():
            flatten_arrays(item, prefix + ('__' if prefix else '') + str(key), result)
    elif isinstance(value, np.ndarray):
        result[prefix] = value
    elif isinstance(value, list) and value and isinstance(value[0], (int, float)):
        result[prefix] = np.asarray(value)
    return result


def json_ready(value):
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_ready(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='A new or empty directory; existing results are never overwritten')
    parser.add_argument('--track-pt', type=float, default=2.)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error('Output directory is not empty')
    data, provenance = load_inputs(args.input_dir, args.track_pt)
    checks = validate_metrics(data)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    from responseDiagnosticsDrawing import render
    plots = render(data, args.output_dir)
    np.savez_compressed(args.output_dir / 'response_diagnostics.npz', **flatten_arrays(data))
    with (args.output_dir / 'migration_metrics.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['fold', 'channel', 'jet_pt_low', 'jet_pt_high',
                         'matched_purity', 'matched_stability',
                         'feed_in_low', 'feed_in_analysis', 'feed_in_high',
                         'feed_out_low', 'feed_out_analysis', 'feed_out_high'])
        for fold in data['folds']:
            for channel in ['jet'] + data['powers']:
                metric = data['metrics'][fold][channel]
                for position, index in enumerate(data['reported_bins']):
                    writer.writerow([fold, channel, data['jet_edges'][index], data['jet_edges'][index + 1],
                                     metric['purity'][index], metric['stability'][index],
                                     *metric['feed_in'][position], *metric['feed_out'][position]])
    notes = [
        'Inclusive OO simulation, R=0.4; particle pT threshold 2 GeV; weight exponents N=1 and N=2.',
        'These EEC responses use generator-particle pairs and encode jet-pT migration, not reconstructed-track angular resolution.',
        'ROOT matrix X=reco and Y=truth; NPZ arrays are [reco, truth]. Matrix plots display truth on x and reco on y.',
        'Stored response contents are weighted bin sums, without bin-width normalization. Counts labels mean weighted counts/sums.',
        'Conditional responses divide each truth column by its matched sum over ALL stored reco bins.',
        'Matched purity = diagonal / matched row sum. Matched stability = diagonal / matched column sum.',
        'EEC jet-bin purity, stability and feed fractions first sum over both angular dimensions.',
        'Feed-in: truth groups <40, 40-160, >=160 GeV within each selected reco bin. Feed-out: reco groups within each selected truth bin.',
        'Feed fractions are conditional on matched pairs/jets. Reco jets below 40 GeV are excluded at production; zero low-reco feed-out is not zero inefficiency.',
        'Truth guard bins 5-40 and >=160 GeV are retained in all denominators. The top stored bin extends to 5020 GeV.',
        'Angular shape panels divide each jet block by its own sum and by DeltaR bin width; they compare unit-area shapes, not absolute EEC yields.',
        'Measured and truth spectra can include unmatched contributions and need not equal response marginals.',
        'Correlation Cij/sqrt(Cii*Cjj) uses stored within-jet second moments; it is not a full event covariance or a mean-subtracted Pearson estimate.',
        'Zero-support ratios/conditional columns are NaN in NPZ, null in JSON, and masked in figures. No interpolation, smoothing, or zero-denominator replacement.',
        'Training/validation comparisons and ratios are descriptive; no statistical closure claim or uncertainty bars without the appropriate covariance.',
        'No unfolding is performed. All four source files were checked against the finalization receipt; histogram underflow/overflow contents were checked to be zero.',
    ]
    (args.output_dir / 'README.txt').write_text('\n'.join(notes) + '\n\nReproduce with:\n'
        'python3 plotOOResponseDiagnostics.py --input-dir .. --output-dir ../response_diagnostics_rerun\n')
    scripts = {}
    for name in ['plotOOResponseDiagnostics.py', 'responseDiagnosticsDrawing.py']:
        path = Path(__file__).parent / name
        shutil.copy2(path, args.output_dir / name)
        scripts[name] = sha256(path)
    summary = {'created_utc': datetime.now(timezone.utc).isoformat(),
               'input_dir': str(args.input_dir.absolute()),
               'selection': {'system': 'OO', 'centrality': 'inclusive', 'jet_radius': 0.4,
                             'particle_pt_min_gev': args.track_pt, 'weight_powers': data['powers'],
                             'jet_edges': data['jet_edges'], 'delta_r_edges': data['dr_edges'],
                             'reported_bins': data['reported_bins']},
               'conventions': notes, 'inputs': provenance, 'validation': checks,
               'scripts_sha256': scripts, 'plots': plots,
               'metrics': {fold: {channel: {key: value for key, value in metrics.items()
                            if key in ['purity', 'stability', 'feed_in', 'feed_out']}
                            for channel, metrics in data['metrics'][fold].items()}
                            for fold in data['folds']}}
    summary['outputs'] = {path.name: {'bytes': path.stat().st_size, 'sha256': sha256(path)}
                          for path in sorted(args.output_dir.iterdir()) if path.is_file()}
    (args.output_dir / 'plot_manifest.json').write_text(json.dumps(json_ready(summary), indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': 'COMPLETE', 'plots': len(plots),
                      'output_dir': str(args.output_dir), 'validation': checks}, indent=2))


if __name__ == '__main__':
    main()
