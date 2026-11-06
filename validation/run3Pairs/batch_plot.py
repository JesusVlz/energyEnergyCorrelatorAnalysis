#!/usr/bin/env python3
"""Use the established pair plotter with reduced physical-event statistics.

This adapter does not reconstruct event vectors or change the correction/plotting
procedure. Batch reductions must supply the same second moments and normalized
shape influence moments that plot_closure.py obtains from physical-event vectors.
All statistical errors remain conditional on the fixed training correction map.
"""
import argparse
from collections import defaultdict
import copy
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import plot_closure as plots


def load_json(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return json.load(handle)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def statistics_index(rows, spectra, label):
    """Require measured moments; only empty zero-yield rows may default to zero."""
    required = [field + suffix for field in plots.FIELDS for suffix in ('2', '_gen')]
    index = {}
    for row in rows:
        key = tuple(row['key'])
        if len(key) != 6 or key in index:
            raise ValueError('{} has an invalid or duplicate key {}'.format(label, key))
        values = {}
        for name in required:
            value = row.get(name)
            if not plots.finite(value):
                raise ValueError('{} missing/nonfinite {} for {}'.format(label, name, key))
            if name.endswith('2') and value < 0:
                raise ValueError('{} negative second moment {} for {}'.format(label, name, key))
            values[name] = value
        if not math.isclose(values['gen_gen'], values['gen2'], rel_tol=1e-10, abs_tol=1e-24):
            raise ValueError('{} truth second moments disagree for {}'.format(label, key))
        index[key] = values
    for row in spectra:
        key = tuple(row['key'])
        if key in index:
            continue
        yields = [row.get(plots.EVENT_FIELDS.get(field, field), row.get(field)) for field in plots.FIELDS]
        if row.get('status') != 'EMPTY' or any(value not in (None, 0) for value in yields):
            raise ValueError('{} missing moments for populated bin {}'.format(label, key))
        index[key] = {name: 0. for name in required}
    return index


def prepare_statistics(result):
    if result.get('holdout_event_vectors'):
        raise ValueError('Batch plotting expects reduced moments, not retained event vectors')
    if 'event_statistics' not in result or 'shape_statistics' not in result:
        raise ValueError('Batch closure requires event_statistics and shape_statistics')
    counts = result.get('event_counts')
    if not isinstance(counts, dict) or not all(name in counts for name in ('selected', 'train', 'holdout', 'raw')):
        raise ValueError('Batch closure requires actual raw/selected/train/holdout event_counts')
    shape_result, unused_empty_vector_statistics, shape_support = plots.normalized_shapes(result)
    absolute = statistics_index(result['event_statistics'], result['holdout_spectra'], 'event_statistics')
    shape = statistics_index(result['shape_statistics'], shape_result['holdout_spectra'], 'shape_statistics')
    return absolute, shape_result, shape, shape_support


def bounded_resolution_inputs(paths, per_file, directory):
    """Materialize only a bounded representative sample in worker-local scratch."""
    outputs, provenance = [], []
    for index, path in enumerate(paths):
        output = Path(directory) / ('resolution_events_{}.jsonl'.format(index))
        opener = gzip.open if str(path).endswith('.gz') else open
        rows = 0
        with opener(path, 'rt') as source, output.open('w') as destination:
            for line in source:
                if not line.strip():
                    continue
                destination.write(line if line.endswith('\n') else line + '\n')
                rows += 1
                if rows >= per_file:
                    break
        outputs.append(output)
        provenance.append({'path': str(path), 'events_copied': rows, 'limit_per_file': per_file})
    return outputs, provenance


def render(ROOT, result, output, input_provenance, selected_jet=None,
           tables_only=False, events=None, resolution_events_per_file=10000):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    result = copy.copy(result)
    systems = {row['key'][0] for row in result['map'] + result['jet_independent_map']}
    result['holdout_spectra'] = [row for row in result.get('holdout_spectra', []) if row['key'][0] in systems]
    report = {
        **input_provenance,
        'adapter_sha256': file_sha256(Path(__file__)),
        'core_plotter_sha256': file_sha256(Path(plots.__file__)),
        'ROOT_version': ROOT.gROOT.GetVersion(),
        'event_counts': result.get('event_counts'),
        'statistics_source': 'Reduced sums of physical-event products and normalized-shape influences; no pair-Poisson errors',
        'uncertainties': 'Physical-event statistical errors conditional on fixed training map and common jet normalization; finite-training-map uncertainty remains pending',
        'shape_uncertainties': 'Joint physical-event delta-method covariance includes each curve normalization and numerator/truth covariance; conditional on fixed training map',
        'open_marker_policy': 'Open squares show supported pair contributions only; incomplete pair shapes are omitted, never partially renormalized',
        'AN_reference': result.get('AN_reference'),
        'correction_version': result.get('correction_version'),
        'shape_convention': 'Established AN interval [0.008,0.4); full stored angular support retained in ROOT',
        'tables': plots.write_tables(ROOT, result, output),
    }
    if not tables_only:
        absolute, shape_result, shape_statistics, shape_support = prepare_statistics(result)
        rootfile = ROOT.TFile(str(output / 'closure.root'), 'RECREATE')
        try:
            report['closure'] = plots.draw_closure(ROOT, result, output, rootfile, selected_jet,
                                                  statistics_override=absolute)
            report['shape_closure'] = plots.draw_closure(ROOT, shape_result, output, rootfile, selected_jet,
                                                        shape=True, statistics_override=shape_statistics)
            report['shape_support'] = shape_support
            if events:
                with tempfile.TemporaryDirectory(prefix='eec_pair_resolution_',
                                                 dir=os.environ.get('_CONDOR_SCRATCH_DIR')) as scratch:
                    paths, provenance = bounded_resolution_inputs(events, resolution_events_per_file, scratch)
                    report['resolution'] = plots.resolution(ROOT, result, paths, output, rootfile)
                    report['resolution_inputs'] = provenance
                report['resolution_method'] = 'Existing greedy charge/dR/pT matching on physical common-reco-jet tracks, 40-80 GeV, held-out events; single-track factor once; no new fit'
            rootfile.cd()
            plots.put_json(ROOT, 'Provenance', report)
        finally:
            rootfile.Close()
    (output / 'support_and_closure.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='Merged batch validation JSON, optionally gzip compressed')
    parser.add_argument('output', type=Path)
    parser.add_argument('--jet-bin', type=int, choices=(0, 1, 2))
    parser.add_argument('--first-then-all', action='store_true', help='Render first_40_80/ before all/, using the same global correction')
    parser.add_argument('--tables-only', action='store_true')
    parser.add_argument('--events', type=Path, nargs='*', default=[], help='Optional representative holdout JSONL or gzip caches for existing resolution plots')
    parser.add_argument('--resolution-events-per-file', type=int, default=10000)
    args = parser.parse_args()
    if args.first_then_all and (args.jet_bin is not None or args.tables_only):
        parser.error('--first-then-all cannot be combined with --jet-bin or --tables-only')
    if args.resolution_events_per_file < 1:
        parser.error('--resolution-events-per-file must be positive')
    result = load_json(args.input)
    provenance = {'input': str(args.input.resolve()), 'input_sha256': file_sha256(args.input)}
    ROOT = plots.root_setup()
    outputs = [(args.output / 'first_40_80', 0), (args.output / 'all', None)] if args.first_then_all else [(args.output, args.jet_bin)]
    summaries = []
    for output, selected_jet in outputs:
        report = render(ROOT, result, output, provenance, selected_jet, args.tables_only,
                        args.events, args.resolution_events_per_file)
        summaries.append({'output': str(output.resolve()), 'tables': len(report['tables']),
                          'closure_plots': len(report.get('closure', [])),
                          'shape_plots': len(report.get('shape_closure', [])),
                          'resolution_plots': len(report.get('resolution', []))})
    print(json.dumps({'status': 'COMPLETE_CONDITIONAL_STATISTICS', 'products': summaries}))


if __name__ == '__main__':
    main()
