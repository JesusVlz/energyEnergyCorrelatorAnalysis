#!/usr/bin/env python3
"""Make Python validation plots for OO GenGen backgrounds on the stored EEC scale.

The inclusive input is used without a centrality choice, like pp and pPb.
Coefficients match EECHistogramManager::SubtractBackground(0/1). Projected S1
already merges the two mixed origins. Shared-event covariance is unavailable.
Defaults: finite jet bins inside 40-200 GeV and charged-particle pT > 2 GeV.
For C++ closure figures use compareDifferentEECpairingsRun3.C; see README_Run3Closure.md.
"""
import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

# Cache floating-point limits before ROOT changes the process FP environment.
np.finfo(np.float64)

PAIR_NAMES = dict(enumerate(('SignalConePair', 'SignalReflectedConePair', 'ReflectedConePair',
    'SignalMixedConePair', 'MixedConePair', 'SignalSecondMixedConePair', 'MixedMixedConePair',
    'SecondMixedConePair', 'SignalPerpendicularConePair', 'PerpendicularConePair',
    'SignalSecondPerpendicularConePair', 'PerpendicularPerpendicularConePair', 'SecondPerpendicularConePair')))


def coefficients():
    result = {name: {(0, origin): 1.} for name, origin in (('T', 4), ('PP', 0), ('PH', 1), ('HH', 3))}
    result['Btrue'] = {(0, 1): 1., (0, 3): 1.}
    for prefix, first, second, same1, same2, cross in (('mix', 3, 5, 4, 7, 6), ('perp', 8, 10, 9, 12, 11)):
        for name, pairing in (('J1', first), ('J2', second), ('same1', same1), ('same2', same2), ('cross', cross)):
            result[prefix + '_' + name] = {(pairing, 4): 1.}
        result[prefix + '_jet_avg'] = {(first, 4): .5, (second, 4): .5}
        result[prefix + '_same_avg'] = {(same1, 4): .5, (same2, 4): .5}
        result[prefix + '_PMavg'] = {(first, 1): .5, (second, 1): .5}
        result[prefix + '_HMavg_minus_cross'] = {(first, 3): .5, (second, 3): .5, (cross, 4): -1.}
        background = {(first, 4): .5, (second, 4): .5, (same1, 4): .5, (same2, 4): .5, (cross, 4): -1.}
        result['B' + prefix] = background
        result['S' + prefix] = {(0, 4): 1., **{key: -weight for key, weight in background.items()}}
    return result


COEFFICIENTS = coefficients()


def combine(raw):
    values = {name: sum(w * raw[key][0] for key, w in terms.items()) for name, terms in COEFFICIENTS.items()}
    variances = {name: sum(w*w * raw[key][1] for key, w in terms.items()) for name, terms in COEFFICIENTS.items()}
    return values, variances


def read_card(source):
    directory = source.Get('JCard')
    if not directory:
        raise ValueError('Missing JCard')
    card = {}
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom('TVectorT<float>'):
            card[key.GetName()] = [float(obj.GetMatrixArray()[i]) for i in range(obj.GetNrows())]
        elif obj.InheritsFrom('TObjString'):
            card[key.GetName()] = str(obj.GetString())
    return card


def number(value):
    return format(value, 'g').replace('-', 'm').replace('.', 'p')


def read_histogram(source, pairing, origin, cent, jet, track, njet):
    token = PAIR_NAMES[pairing]
    path = 'energyEnergyCorrelator/{0}/energyEnergyCorrelator{0}_C{1}T{2}{3}{4}'.format(
        token, cent, track, 'J{}'.format(jet) if jet < njet else '', 'S{}'.format(origin) if origin < 4 else '')
    hist = source.Get(path)
    if not hist or not hist.InheritsFrom('TH1'):
        raise ValueError('Missing histogram: ' + path)
    count = hist.GetNbinsX()
    edges = np.array([hist.GetBinLowEdge(i) for i in range(1, count + 2)])
    values = np.array([hist.GetBinContent(i) for i in range(1, count + 1)])
    variances = np.array([hist.GetBinError(i)**2 for i in range(1, count + 1)])
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(variances)):
        raise ValueError('Nonfinite content/error: ' + path)
    return edges, values, variances


def load_channels(ROOT, filename, scope):
    source = ROOT.TFile.Open(str(filename), 'READ')
    if not source or source.IsZombie():
        raise ValueError('Cannot open ' + str(filename))
    card = read_card(source)
    if card['DataType'] != [10.] or card['McCorrelationType'] != [3.]:
        raise ValueError('Requires OO GenGen MC')
    if int(card['BackgroundMethods'][0]) & 6 != 6 or len(card['WeightExponent']) != 1:
        raise ValueError('Requires mixed/perpendicular categories and a single projected exponent')
    centralities, jets, tracks = (card[key] for key in ('CentralityBinEdges', 'JetPtBinEdgesEEC', 'TrackPtBinEdgesEEC'))
    required = sorted(set(key for terms in COEFFICIENTS.values() for key in terms))
    # Inclusive uses all native slots, including sentinel/peripheral bins: no centrality selection.
    selections = [list(range(len(centralities)-1))] if scope == 'integrated' else [[i] for i in range(len(centralities)-1)]
    channels = []
    for selection in selections:
        cent = selection[0]
        low, high = centralities[cent], centralities[selection[-1]+1]
        if scope == 'integrated':
            label, cname = 'Inclusive OO', 'inclusive'
        elif high <= 0:
            label, cname = 'Sentinel/unassigned centrality (-1 to 0)', 'Csentinel'
        else:
            label = 'Centrality {:g}\N{EN DASH}{:g}%'.format(low, high)
            cname = 'C{}-{}'.format(number(low), number(high))
        for jet in range(len(jets)):
            inclusive = jet == len(jets) - 1
            borders = [card['MinJetPtCut'][0], card['MaxJetPtCut'][0]] if inclusive else jets[jet:jet+2]
            jname = 'Jinclusive' if inclusive else 'J{}-{}'.format(*map(number, borders))
            jlabel = 'Inclusive selected jets (no jet-axis restriction)' if inclusive else r'{:g} < jet $p_T$ < {:g} GeV'.format(*borders)
            for track, threshold in enumerate(tracks[:-1]):
                raw, template = {}, None
                for pairing, origin in required:
                    accumulated, accumulated_var = None, None
                    for native_cent in selection:
                        edges, values, variances = read_histogram(source, pairing, origin, native_cent, jet, track, len(jets)-1)
                        if template is not None and not np.array_equal(template, edges):
                            raise ValueError('Histogram axes differ within a channel')
                        template = edges
                        accumulated = values.copy() if accumulated is None else accumulated + values
                        accumulated_var = variances.copy() if accumulated_var is None else accumulated_var + variances
                    raw[pairing, origin] = accumulated, accumulated_var
                values, variances = combine(raw)
                residual = values['T'] - values['PP'] - values['PH'] - values['HH']
                relative = float(np.max(np.abs(residual)) / max(float(np.max(np.abs(values['T']))), 1e-300))
                if relative > 1e-8:
                    raise ValueError('Subevent sum failed: {} C{} J{} T{}: {}'.format(filename, selection, jet, track, relative))
                identifier = '{}_{}_T{}_n{}'.format(cname, jname, number(threshold), number(card['WeightExponent'][0]))
                channels.append(dict(id=identifier, scope=scope, centrality=[low, high], centrality_indices=selection,
                    centrality_index=cent, centrality_name=cname, label=label, jet=borders, jet_label=jlabel, jet_index=jet,
                    jet_inclusive=inclusive, track_threshold=threshold, exponent=card['WeightExponent'][0],
                    radius=card['JetRadius'][0], edges=edges, values=values, variances=variances,
                    status='populated' if np.any(values['T'] != 0) else 'empty',
                    decomposition_max_relative_residual=relative))
    source.Close()
    return card, channels


def channel_summary(channel):
    result = {k: v for k, v in channel.items() if k not in ('edges', 'values', 'variances')}
    edges, values = channel['edges'], channel['values']
    result['integrals'] = {}
    for view, low, high in (('legacy', .008, .39), ('full', 0., 2*channel['radius'])):
        mask = (edges[:-1] >= low-1e-8) & (edges[1:] <= high+1e-8)
        bins = np.flatnonzero(mask)
        integrals = {name: float(np.dot(y[mask], np.diff(edges)[mask])) for name, y in values.items()}
        ratios = {estimate + '_over_' + truth: integrals[estimate]/integrals[truth] if integrals[truth] > 0 else None
                  for estimate, truth in (('Bmix', 'Btrue'), ('Bperp', 'Btrue'), ('Smix', 'PP'), ('Sperp', 'PP'))}
        result['integrals'][view] = dict(requested_range=[low, high],
            stored_bin_range=[float(edges[bins[0]]), float(edges[bins[-1]+1])] if len(bins) else None,
            values=integrals, ratios=ratios)
    result['empty_method_terms'] = {prefix: [term for term in ('J1', 'J2', 'same1', 'same2', 'cross')
        if not np.any(values[prefix + '_' + term] != 0)] for prefix in ('mix', 'perp')}
    result['nonpositive_signal_bins'] = {name: int(np.count_nonzero(values[name] <= 0)) for name in ('PP', 'Smix', 'Sperp')}
    return result


def focus_channels(channels, jet_range=(40., 200.), track_threshold=2.):
    """Select finite native jet bins and the exact requested track threshold."""
    low, high = jet_range
    if low >= high:
        raise ValueError('Jet range must have increasing boundaries')
    return [c for c in channels if not c['jet_inclusive']
            and c['jet'][0] >= low and c['jet'][1] <= high
            and abs(c['track_threshold'] - track_threshold) < 1e-8]


def draw_job(job):
    from ooBackgroundDrawing import draw_channel
    channel, destination = job
    draw_channel(channel, destination)
    return channel['id']


def concatenate(paths, destination):
    if paths:
        subprocess.run(['pdfunite', *map(str, paths), str(destination)], check=True)


def write_root(ROOT, channels, destination):
    output = ROOT.TFile.Open(str(destination), 'RECREATE')
    for channel in channels:
        output.mkdir(channel['id']).cd()
        for name, values in channel['values'].items():
            histogram = ROOT.TH1D(name, name, len(values), array('d', channel['edges']))
            histogram.SetDirectory(0)
            for index, value in enumerate(values, 1):
                histogram.SetBinContent(index, value)
                histogram.SetBinError(index, math.sqrt(channel['variances'][name][index-1]))
            histogram.Write()
        output.cd()
    ROOT.TNamed('uncertainties', 'Saved histogram diagonal errors; shared-event covariance unavailable').Write()
    ROOT.TNamed('normalization', 'Stored EEC density; no independent rescaling; regular bins only').Write()
    output.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--integrated', type=Path, required=True)
    parser.add_argument('--split', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--jet-pt-range', type=float, nargs=2, default=(40., 200.), metavar=('MIN', 'MAX'))
    parser.add_argument('--track-threshold', type=float, default=2.)
    parser.add_argument('--resume', action='store_true', help='Skip complete plots in an existing output directory')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=args.resume)
    import ROOT
    ROOT.gROOT.SetBatch(True)
    cards, channels, inputs = {}, [], []
    for scope, filename in (('integrated', args.integrated), ('split', args.split)):
        print('Reading', scope, filename, flush=True)
        card, current = load_channels(ROOT, filename, scope)
        cards[scope] = card
        channels.extend(current)
        inputs.append(dict(scope=scope, path=str(filename.resolve()), bytes=filename.stat().st_size,
                           sha256=hashlib.sha256(filename.read_bytes()).hexdigest()))
    channels = focus_channels(channels, args.jet_pt_range, args.track_threshold)
    if not channels:
        raise ValueError('No stored channels match the requested jet and track selection')
    selection = dict(jet_pt_range=list(args.jet_pt_range), track_threshold=args.track_threshold, inclusive_jet_slot=False)
    manifest_path = output/'manifest.json'
    previous = json.loads(manifest_path.read_text()) if args.resume and manifest_path.exists() else None
    if previous is not None and previous['inputs'] != inputs:
        raise ValueError('Cannot resume: input hashes differ')
    if previous is not None and previous.get('selection') != selection:
        raise ValueError('Cannot resume: plotting selection differs')
    populated = [c for c in channels if c['status'] == 'populated']
    print('Read {} channels; {} populated'.format(len(channels), len(populated)), flush=True)
    write_root(ROOT, channels, output/'derived_histograms.root')
    summary = dict(cards=cards, channels=[channel_summary(c) for c in channels],
        inclusive_policy='all native centrality slots; no centrality choice or label in inclusive plots',
        normalization='stored EEC density; common scale, no independent rescaling',
        uncertainties='saved histogram diagonal errors; shared-event covariance unavailable', selection=selection, grids=False)
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    with (output/'integrals.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('channel', 'view', 'metric', 'integral', 'ratio_to_truth'))
        for channel in summary['channels']:
            for view, entry in channel['integrals'].items():
                for metric in ('T', 'PP', 'PH', 'HH', 'Btrue', 'Bmix', 'Bperp', 'Smix', 'Sperp'):
                    truth = 'Btrue' if metric.startswith('B') else 'PP'
                    writer.writerow((channel['id'], view, metric, entry['values'][metric],
                                     entry['ratios'].get(metric+'_over_'+truth)))
    snapshot = output/'source_snapshot'
    snapshot.mkdir(exist_ok=True)
    sources = []
    for name in ('plotOOMCBackgroundStudies.py', 'ooBackgroundDrawing.py'):
        original = Path(__file__).with_name(name)
        shutil.copyfile(original, snapshot/name)
        sources.append(dict(path=str(original.resolve()), sha256=hashlib.sha256(original.read_bytes()).hexdigest()))
    manifest = dict(status='rendering', created_utc=datetime.now(timezone.utc).isoformat(), command=sys.argv,
        ROOT=ROOT.gROOT.GetVersion(), inputs=inputs, sources=sources, selection=selection, grids=False, channel_count=len(channels),
        populated_channels=len(populated), empty_channels=len(channels)-len(populated), input_files_modified=False,
        maximum_decomposition_relative_residual=max(c['decomposition_max_relative_residual'] for c in channels))
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    reuse_plots = args.resume and previous is not None and previous.get('sources') == sources
    core = []
    for channel in populated:
        if channel['track_threshold'] == 2. and channel['jet_index'] in (0, 1, 2, 3):
            core.append(dict(channel, id=channel['id']+'_core', display_range=[.008, .39]))
    jobs = []
    for channel in core + sorted(populated, key=lambda c: (c['track_threshold'] != 2., c['jet_inclusive'], c['jet_index'])):
        destination = output/channel['scope']/channel['id']
        destination.parent.mkdir(exist_ok=True)
        if reuse_plots and all(Path(str(destination)+suffix).exists() for suffix in ('.pdf', '.png', '_components.png')):
            continue
        jobs.append((channel, destination))
    with ProcessPoolExecutor(max_workers=max(1, args.workers), mp_context=multiprocessing.get_context('spawn')) as executor:
        for count, identifier in enumerate(executor.map(draw_job, jobs), 1):
            if count % 8 == 0 or count == len(jobs):
                print('Rendered {}/{}: {}'.format(count, len(jobs), identifier), flush=True)
    from ooBackgroundDrawing import draw_overlays
    overlays = output/'comparisons'
    overlays.mkdir(exist_ok=True)
    comparison_pdfs = []
    for group, keys in (('jet', ('scope', 'centrality_index', 'track_threshold')),
                        ('track', ('scope', 'centrality_index', 'jet_index')),
                        ('centrality', ('jet_index', 'track_threshold'))):
        groups = {}
        for channel in populated:
            if group == 'jet' and channel['jet_inclusive']:
                continue
            if group == 'centrality' and channel['scope'] == 'split' and channel['centrality'][1] <= 0:
                continue
            groups.setdefault(tuple(channel[key] for key in keys), []).append(channel)
        for index, current in enumerate(groups.values()):
            if len(current) < 2:
                continue
            destination = overlays/(group+'_'+current[0]['scope']+'_'+str(index))
            if not (reuse_plots and destination.with_suffix('.pdf').exists() and destination.with_suffix('.png').exists()):
                draw_overlays(current, destination, group)
            comparison_pdfs.append(destination.with_suffix('.pdf'))
        print('Rendered', group, 'comparisons', flush=True)
    for scope in cards:
        concatenate([output/scope/(c['id']+'.pdf') for c in populated if c['scope'] == scope],
                    output/(scope+'_background_atlas.pdf'))
    concatenate([output/c['scope']/(c['id']+'.pdf') for c in core], output/'background_overview.pdf')
    concatenate(comparison_pdfs, output/'comparison_atlas.pdf')
    lines = ['OO GenGen background studies', '', 'Start with background_overview.pdf, then integrated_background_atlas.pdf,',
        'split_background_atlas.pdf and comparison_atlas.pdf. Individual PDFs have two',
        'pages (closure, then components); PNGs provide previews.',
        'The overview uses the standard 0.008 < Delta r < 0.39 range at pTch>2 GeV',
        'for the four jet intervals between 40 and 200 GeV. Full atlases extend to 2R.', '',
        'Inclusive OO follows the pp/pPb plotting strategy: no centrality choice and no',
        'centrality bounds in inclusive figures. All native slots, including sentinel and',
        '90-100% slots when present, contribute to the inclusive sample.',
        'Split-centrality plots retain their native centrality labels.', '',
        'Coverage: native centrality splits and inclusive OO; finite jet bins inside',
        '{}-{} GeV and pTch>{} GeV. All plots have no grid.'.format(*args.jet_pt_range, args.track_threshold),
        'Only the stored exponent n=1 is available.',
        'Zero-content selections are listed below and retained in derived_histograms.root.',
        'No reflected-cone histograms exist. Mixed and perpendicular components are shown.',
        'Empty stored background terms are annotated; matching/exposure ledgers are absent.',
        'The no-J integrated jet slot is excluded from this focused plot set.',
        'Plots use the supplied GenGen files as-is; no reconstructed-jet correction is',
        'reapplied, and projected cards do not certify the specific Run3 JEC payload.',
        'Signed values are retained; ratios to nonpositive truth are masked. Error bars',
        'use saved histogram diagonal errors, without shared-event covariance.', '',
        'B = (JC1+JC2)/2 - C1C2 + (C1C1+C2C2)/2; extracted signal = T-B.',
        'Truth signal=S0; truth background=S1+S3. Projected S1 merges both mixed origins.',
        'All curves retain the stored scale without independent normalization.',
        'Regular bins only; full view extends to 2R. Summary integrals use complete',
        'stored bins inside each interval and record their actual boundaries.', '',
        'Numerical evidence: summary.json, integrals.csv, derived_histograms.root.',
        'Input provenance and reproduction sources: manifest.json, source_snapshot/.', '', 'Coverage:']
    for channel in channels:
        lines.append('{}: {}{}'.format(channel['id'], channel['status'],
            ' -> '+channel['scope']+'/'+channel['id']+'.pdf' if channel['status']=='populated' else ' (no entries)'))
    (output/'INDEX.txt').write_text('\n'.join(lines)+'\n')
    for original in inputs:
        actual = hashlib.sha256(Path(original['path']).read_bytes()).hexdigest()
        if actual != original['sha256']:
            raise ValueError('Input changed during rendering: '+original['path'])
    manifest.update(status='complete', comparison_plots=len(comparison_pdfs), core_plots=len(core))
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    print('Finished:', output, flush=True)


if __name__ == '__main__':
    main()
