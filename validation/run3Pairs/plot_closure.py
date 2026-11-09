#!/usr/bin/env python3
"""Write AN-style Run-3 pair tables and holdout EEC closure ROOT/PDF products.

The existing JDrawer supplies the split-pad appearance. Statistical bars use
physical-event vectors, conditional on the fixed training map and jet denominator;
no independent-pair Poisson uncertainty or finite-training uncertainty is implied.
Unsupported table cells are NaN. Missing closure points remain absent; open points
explicitly show only the supported reconstructed-pair contribution.
"""
import argparse
from array import array
from collections import defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path
from preparation_contract import write_contract

FIELDS = ('gen', 'raw', 'single', 'pair', 'pair_jet_independent')
EVENT_FIELDS = {'pair': 'pair_supported', 'pair_jet_independent': 'pair_jet_independent_supported'}
LABELS = {'gen': 'Truth, common reco jets', 'raw': 'Raw reco', 'single': 'Single-track corrected',
          'pair': 'Pair corrected', 'pair_jet_independent': 'Jet-independent pair corrected'}
COLORS = {'gen': 1, 'raw': 632, 'single': 600, 'pair': 416 + 2, 'pair_jet_independent': 800 + 7}


def root_setup():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)
    ROOT.gInterpreter.Declare('#include "' + str(Path(__file__).resolve().parents[2] / 'plotting/JDrawer.h') + '"')
    return ROOT


def finite(value):
    return value is not None and math.isfinite(value)


def slug(key):
    return '_'.join(str(x) for x in key)


def put_json(ROOT, name, data):
    ROOT.TObjString(json.dumps(data, sort_keys=True, allow_nan=False)).Write(name)


def write_tables(ROOT, result, output):
    """Legacy names store epsilon_norm, so the runtime applies its reciprocal."""
    pt, jet, dr = (result['binning'][name] for name in ('pt', 'jet', 'dr_support'))
    rows = result['map']
    independent = result['jet_independent_map']
    systems = sorted({row['key'][0] for row in rows + independent})
    outputs = []
    for system in systems:
        filename = output / ('trackPairEfficiencyCorrection_' + system + '_' + result['correction_version'] + '.root')
        rootfile = ROOT.TFile(str(filename), 'RECREATE')
        ROOT.TNamed('Run3PairCorrectionVersion', result['correction_version']).Write()
        write_contract(ROOT, result, rootfile, system)
        put_json(ROOT, 'Provenance', {key: result.get(key) for key in (
            'method', 'control', 'fold_algorithm', 'min_raw_pairs_per_cell_and_control',
            'input_path', 'input_sha256', 'input_files', 'tool_sha256', 'status', 'invalid_cell_policy', 'correction_version', 'AN_reference', 'jet_independent_population', 'jet_independent_trigger_abs_eta_max_exclusive')})
        card = rootfile.mkdir('JCard')
        card.cd()
        vectors = {'DataType': [10 if system == 'OO' else 12],
                   'CentralityBinEdges': [0., 30., 50., 90., 100.] if system == 'OO' else [-1., 100.],
                   'TrackPairPtBinEdges': pt, 'JetPtBinEdgesEEC': jet,
                   'TriggerEtaCut': [result.get('jet_independent_trigger_abs_eta_max_exclusive', 1.6)],
                   'JetEtaCut': [1.6], 'TrackEtaCut': [2.4]}
        for name, values in vectors.items():
            vector = ROOT.TVectorT('float')(1, len(values))
            for i, value in enumerate(values, 1):
                vector[i] = value
            vector.Write(name)
        rootfile.cd()
        counts = defaultdict(int)
        for collection, nearjet in ((independent, False), (rows, True)):
            lookup = {tuple(row['key']): row for row in collection if row['key'][0] == system}
            activity_count = 4 if system == 'OO' else 1
            for activity in range(activity_count):
                for high in range(len(pt) - 1):
                    for low in range(high + 1):
                        for j in (range(len(jet) - 1) if nearjet else (None,)):
                            name = ('trackPairEfficiencyCorrectionCloseToJet' if nearjet else 'trackPairEfficiencyCorrection') + f'_C{activity}T{high}A{low}'
                            if nearjet:
                                name += f'J{j}'
                            hist = ROOT.TH1D(name, ';#DeltaR;normalized single-track-corrected reco / gen', len(dr) - 1, array('d', dr))
                            mask = ROOT.TH1I(name + '_valid', ';#DeltaR;finite supported correction', len(dr) - 1, array('d', dr))
                            ngen = ROOT.TH1D(name + '_genPairs', ';#DeltaR;raw generator pairs', len(dr) - 1, array('d', dr))
                            nreco = ROOT.TH1D(name + '_recoPairs', ';#DeltaR;raw reconstructed pairs', len(dr) - 1, array('d', dr))
                            alpha = ROOT.TH1D(name + '_controlSingleOverGen', ';#DeltaR;control single / gen', len(dr) - 1, array('d', dr))
                            for ibin in range(1, len(dr)):
                                key = (system, 'R04', activity, j, high, low, ibin - 1) if nearjet else (system, 'R04', activity, high, low, ibin - 1)
                                cell = lookup.get(key, {})
                                correction = cell.get('correction')
                                valid = bool(cell.get('valid')) and finite(correction) and correction > 0
                                hist.SetBinContent(ibin, 1. / correction if valid else float('nan'))
                                # Map-replica uncertainty is not evaluated here. Zero
                                # bin errors are bookkeeping, not a precision claim.
                                hist.SetBinError(ibin, 0.)
                                mask.SetBinContent(ibin, int(valid))
                                ngen.SetBinContent(ibin, cell.get('gen_pairs', 0))
                                nreco.SetBinContent(ibin, cell.get('reco_pairs', 0))
                                value = cell.get('large_angle_single_over_gen')
                                alpha.SetBinContent(ibin, value if finite(value) else float('nan'))
                                counts[('jet_dependent' if nearjet else 'jet_independent') + ('_supported' if valid else '_unsupported')] += 1
                            # Flow bins must not silently become unity or zero.
                            hist.SetBinContent(0, float('nan'))
                            hist.SetBinContent(len(dr), float('nan'))
                            for obj in (hist, mask, ngen, nreco, alpha):
                                obj.Write()
        put_json(ROOT, 'SupportSummary', dict(counts))
        ROOT.TNamed('UncertaintyStatus', 'finite-training-map uncertainty not evaluated; histogram errors are not uncertainty estimates').Write()
        rootfile.Close()
        outputs.append({'path': str(filename.resolve()), 'system': system, **counts})
    return outputs


def event_statistics(result):
    # Ratio uncertainty includes numerator/truth correlation within a physical event.
    stats = defaultdict(lambda: defaultdict(float))
    for event in result.get('holdout_event_vectors', []):
        for row in event['channels']:
            key = tuple(row['key'])
            truth = row.get('gen', 0.)
            for field in FIELDS:
                value = row.get(EVENT_FIELDS.get(field, field), 0.)
                if value is None:
                    continue
                stats[key][field + '2'] += value * value
                stats[key][field + '_gen'] += value * truth
    return stats


def ratio_error(stat, field, numerator, denominator):
    if denominator <= 0:
        return 0.
    ratio = numerator / denominator
    variance = stat[field + '2'] + ratio * ratio * stat['gen2'] - 2 * ratio * stat[field + '_gen']
    return math.sqrt(max(0., variance)) / denominator


def draw_closure(ROOT, result, output, rootfile, selected_jet=None, shape=False, statistics_override=None):
    dr, jet = result['binning']['dr_support'], result['binning']['jet']
    groups = defaultdict(list)
    for row in result['holdout_spectra']:
        if selected_jet is None or row['key'][2] == selected_jet:
            groups[tuple(row['key'][:-1])].append(row)
    statistics = statistics_override if statistics_override is not None else event_statistics(result)
    norm = {tuple(row['key']): row.get('weighted_jets', 0.) for row in result.get('jet_event_sums', [])}
    directory = rootfile.mkdir('closure_shape' if shape else 'closure')
    report = []
    for key, rows in sorted(groups.items()):
        rows = sorted(rows, key=lambda row: row['key'][-1])
        system, radius, j, n, selection = key
        channel = slug(key) + ('_shape' if shape else '')
        directory.cd()
        outdir = directory.mkdir(channel)
        outdir.cd()
        denominator = 0. if shape else norm.get((system, radius, j, selection), 0.)
        scale = 1. / denominator if denominator > 0 else 1.
        hists, graphs, ratios, partials, partial_ratios = {}, {}, {}, {}, {}
        maximum = 0.
        positive = []
        unsupported = []
        residuals = []
        single_residuals = []
        ratio_limits = [2.]
        keep = []
        for field in FIELDS:
            hist = ROOT.TH1D(field, ';#DeltaR;EEC / #DeltaR', len(dr) - 1, array('d', dr))
            graph, ratio = ROOT.TGraphErrors(), ROOT.TGraphErrors()
            partial, partial_ratio = ROOT.TGraphErrors(), ROOT.TGraphErrors()
            graph.SetName(field + '_supported_bins')
            ratio.SetName(field + '_over_truth')
            partial.SetName(field + '_partial_pairs')
            partial_ratio.SetName(field + '_partial_pairs_over_truth')
            for obj in (graph, ratio, partial, partial_ratio):
                obj.SetLineColor(COLORS[field]); obj.SetMarkerColor(COLORS[field]); obj.SetLineWidth(2)
                obj.SetMarkerStyle(20 if field != 'gen' else 24); obj.SetMarkerSize(.72)
            partial.SetMarkerStyle(25); partial_ratio.SetMarkerStyle(25)
            for row in rows:
                ibin = row['key'][-1]
                low, high = dr[ibin:ibin + 2]
                width = high - low
                value = row.get(field)
                is_partial = field in EVENT_FIELDS and not finite(value)
                if is_partial and shape:
                    hist.SetBinContent(ibin + 1, float('nan'))
                    continue
                if is_partial:
                    value = row.get(EVENT_FIELDS[field], 0.)
                stat = statistics[tuple(row['key'])]
                y = value * scale / width if finite(value) else 0.
                err = math.sqrt(stat[field + '2']) * scale / width
                hist.SetBinContent(ibin + 1, y if finite(row.get(field)) else float('nan'))
                hist.SetBinError(ibin + 1, err if finite(row.get(field)) else 0.)
                if not finite(value) or low < .008 or high > .4:
                    continue
                x = (low + high) / 2
                if value > 0:
                    target = partial if is_partial else graph
                    point = target.GetN(); target.SetPoint(point, x, y); target.SetPointError(point, width / 2, err)
                    maximum = max(maximum, y + err)
                    positive.append(y)
                if row['gen'] > 0 and field != 'gen':
                    target_ratio = partial_ratio if is_partial else ratio
                    point = target_ratio.GetN()
                    r = value / row['gen']
                    rerr = ratio_error(stat, field, value, row['gen'])
                    target_ratio.SetPoint(point, x, r); target_ratio.SetPointError(point, width / 2, rerr)
                    ratio_limits.append(r + rerr)
                    if field == 'single':
                        single_residuals.append({'dr': [low, high], 'ratio': r, 'conditional_event_error': rerr})
                    if not is_partial and field == 'pair':
                        residuals.append({'dr': [low, high], 'ratio': r, 'conditional_event_error': rerr})
            hist.Write(); graph.Write(); ratio.Write(); partial.Write(); partial_ratio.Write()
            hists[field], graphs[field], ratios[field] = hist, graph, ratio
            partials[field], partial_ratios[field] = partial, partial_ratio
        support = ROOT.TH1D('supported_reco_weight_fraction', ';#DeltaR;supported single-track weight / all single-track weight', len(dr) - 1, array('d', dr))
        for row in rows:
            ibin = row['key'][-1]
            fraction = 1. - row.get('masked_single', 0.) / row['single'] if row['single'] > 0 else float('nan')
            support.SetBinContent(ibin + 1, fraction)
            if row.get('masked_pairs', 0) and dr[ibin] < .4:
                unsupported.append({'dr': dr[ibin:ibin + 2], 'masked_pairs': row['masked_pairs'],
                                    'reco_pairs': row.get('reco_pairs'), 'supported_single_weight_fraction': fraction,
                                    'fallback_pairs': row.get('fallback_pairs', 0)})
        support.Write()
        drawer = ROOT.JDrawer()
        ROOT.SetOwnership(drawer, False)  # JDrawer owns its pads; avoid PyROOT double teardown.
        drawer.SetDefaultAppearanceSplitCanvas(); drawer.SetLogX(True); drawer.SetLogY(True)
        drawer.SetBottomMargin(.13); drawer.SetTitleOffsetX(.8)
        drawer.SetRelativeCanvasSize(1.1, .95)
        drawer.CreateSplitCanvas()
        frame = hists['gen'].Clone('frame_' + channel); frame.Reset()
        xminimum = .008
        frame.GetXaxis().SetRangeUser(xminimum, .4)
        frame.SetMinimum(max(min(positive) / 5, 1e-18) if positive else .1)
        frame.SetMaximum(maximum * 30 if maximum > 0 else 10.)
        ylabel = f'(1/#Sigma) d#Sigma^{{({n})}}/d#DeltaR' if shape else (f'(1/N_{{jet}}) d#Sigma^{{({n})}}/d#DeltaR' if denominator > 0 else f'd#Sigma^{{({n})}}/d#DeltaR')
        drawer.DrawHistogramToUpperPad(frame, '#DeltaR', ylabel, '', 'AXIS')
        legend = ROOT.TLegend(.45, .56, .94, .84); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.032)
        label = 'OO selected MB' if selection == 'MB' else ('pp reference' if selection == 'ppRef' else 'OO ' + selection[1:].replace('_', '-') + '%')
        for field in FIELDS:
            graphs[field].Draw('PZ SAME')
            partials[field].Draw('PZ SAME')
            legend.AddEntry(graphs[field], LABELS[field], 'p')
        legend.Draw()
        note = ROOT.TLatex(); note.SetNDC(); note.SetTextFont(42); note.SetTextSize(.033)
        note.DrawLatex(.22, .94, f'{label}, R=0.4, {jet[j]:g} < p_{{T,jet}} < {jet[j+1]:g} GeV, n={n}')
        note.SetTextSize(.028)
        note.DrawLatex(.22, .89, 'p_{T,track}>2 GeV, |#eta|<2.4; independent validation events')
        if shape:
            note.DrawLatex(.22, .49, 'Unit area: 0.008 #leq #DeltaR < 0.4; incomplete pair shapes omitted')
        elif any(graph.GetN() for graph in partials.values()):
            note.DrawLatex(.22, .49, 'Open squares: supported pairs only; incomplete correction')
        drawer.SetLogY(False)
        ratioframe = frame.Clone('ratioframe_' + channel); ratioframe.Reset(); ratioframe.SetMinimum(0.); ratioframe.SetMaximum(max(ratio_limits) * 1.05)
        drawer.DrawHistogramToLowerPad(ratioframe, '#DeltaR', 'Reco / truth', '', 'AXIS')
        line = ROOT.TLine(xminimum, 1., .4, 1.); line.SetLineStyle(2); line.Draw()
        for field in FIELDS[1:]:
            ratios[field].Draw('PZ SAME'); partial_ratios[field].Draw('PZ SAME')
        canvas = ROOT.gPad.GetCanvas(); canvas.SetName(channel)
        filename = output / (channel + '.pdf')
        canvas.SaveAs(str(filename)); canvas.Write('canvas')
        keep.extend((frame, ratioframe, legend, note, line, canvas, drawer))
        report.append({'channel': list(key), 'normalization': 'unit area [0.008,0.4)' if shape else 'common selected weighted jets', 'weighted_jets': denominator,
                       'unsupported_angular_bins': unsupported,
                       'largest_single_track_residuals': sorted(single_residuals, key=lambda x: abs(x['ratio'] - 1), reverse=True)[:3],
                       'largest_supported_pair_residuals': sorted(residuals, key=lambda x: abs(x['ratio'] - 1), reverse=True)[:3],
                       'pdf': str(filename.resolve())})
        # ROOT canvases are streamed before close; retained Python objects survive SaveAs.
        canvas.Close()
    return report



def normalized_shapes(result):
    """Standard AN-style unit-area shapes and joint event influence covariance.

    The literal .008 boundary must exist. A missing pair bin invalidates the
    entire shape; the supported subset is never renormalized to unit area.
    """
    import copy
    dr = result['binning']['dr_support']
    if .008 not in dr or .4 not in dr:
        raise ValueError('Standard shape display requires stored edges .008 and .4; rederive without interpolating counts')
    selected_bins = {i for i in range(len(dr)-1) if .008 <= dr[i] and dr[i+1] <= .4}
    groups = defaultdict(list)
    for row in result['holdout_spectra']:
        if row['key'][-1] in selected_bins:
            groups[tuple(row['key'][:-1])].append(row)
    totals, shapes, diagnostics = {}, {}, []
    for channel, rows in groups.items():
        totals[channel] = {}
        reasons = {}
        for field in FIELDS:
            missing = any(row.get(field) is None and row.get('status') != 'EMPTY' for row in rows)
            total = sum(row.get(field) or 0. for row in rows)
            totals[channel][field] = total if not missing and total > 0 else None
            if missing:
                reasons[field] = 'unsupported populated bin in normalization interval'
            elif total <= 0:
                reasons[field] = 'nonpositive integral in normalization interval'
        for row in rows:
            key = tuple(row['key'])
            shapes[key] = {field: (row.get(field) or 0.) / totals[channel][field] if totals[channel][field] else None for field in FIELDS}
        diagnostics.append({'channel': list(channel), 'integrals': totals[channel], 'suppressed_curves': reasons})
    stats = defaultdict(lambda: defaultdict(float))
    for event in result.get('holdout_event_vectors', []):
        event_groups = defaultdict(dict)
        for row in event['channels']:
            if row['key'][-1] in selected_bins:
                event_groups[tuple(row['key'][:-1])][row['key'][-1]] = row
        for channel, bins in event_groups.items():
            event_totals = {field: sum(row.get(EVENT_FIELDS.get(field, field), 0.) for row in bins.values()) for field in FIELDS}
            for ibin in selected_bins:
                key = (*channel, ibin)
                if key not in shapes:
                    continue
                row = bins.get(ibin, {})
                influence = {}
                for field in FIELDS:
                    total = totals[channel][field]
                    value = row.get(EVENT_FIELDS.get(field, field), 0.)
                    influence[field] = (value - shapes[key][field] * event_totals[field]) / total if total else 0.
                for field in FIELDS:
                    stats[key][field + '2'] += influence[field] ** 2
                    stats[key][field + '_gen'] += influence[field] * influence['gen']
    normalized = {'binning': result['binning'], 'holdout_spectra': [], 'jet_event_sums': []}
    for row in result['holdout_spectra']:
        key = tuple(row['key'])
        if key not in shapes:
            continue
        entry = copy.copy(row)
        for field in FIELDS:
            entry[field] = shapes[key][field]
        # Truth/raw/single absent only for an entirely empty channel. Zero is
        # a legitimate empty histogram; pair curves remain explicitly absent.
        for field in ('gen', 'raw', 'single'):
            if entry[field] is None:
                entry[field] = 0.
        single_total = totals[key[:-1]]['single']
        entry['masked_single'] = row.get('masked_single', 0.) / single_total if single_total else 0.
        entry['pair_supported'] = 0.
        entry['pair_jet_independent_supported'] = 0.
        normalized['holdout_spectra'].append(entry)
    return normalized, stats, diagnostics


def matching_pairs(jet):
    """Existing ConstructParticleResponses greedy charge/dR/pT matching, exactly."""
    reco = sorted(jet['reco'], key=lambda t: (t['pt'], t['eta'], t['phi'], t.get('correction', 1.), t['charge']), reverse=True)
    truth = jet['gen']
    used, matches = set(), []
    for track in reco:
        candidates = []
        for i, particle in enumerate(truth):
            if i in used or track['charge'] != particle['charge']:
                continue
            distance = math.hypot(track['eta'] - particle['eta'], math.remainder(track['phi'] - particle['phi'], 2 * math.pi))
            dpt = abs(track['pt'] - particle['pt'])
            if distance < .05 and dpt < .5 * track['pt']:
                candidates.append((dpt, i))
        if candidates:
            _, i = min(candidates); used.add(i); matches.append((track, truth[i]))
    return itertools.combinations(matches, 2)


def resolution(ROOT, result, event_paths, output, rootfile):
    from pair_counts import fold, separation
    dr = result['binning']['dr_support']
    ptproduct = [4,6,9,12,16,20,25,30,35,40,45,50,60,70,80,90,100,120,140,180,240,400,600,2000,90000]
    histograms = {}
    stats = defaultdict(lambda: defaultdict(float))
    for system in ('OO', 'pp'):
        histograms[system] = {
            'angular': ROOT.TH2D(system + '_angular', ';#DeltaR (reco);#DeltaR (gen)', len(dr)-1, array('d', dr), len(dr)-1, array('d', dr)),
            'momentum': ROOT.TH2D(system + '_momentum', ';p_{T,1}p_{T,2} (reco);p_{T,1}p_{T,2} (gen)', len(ptproduct)-1, array('d', ptproduct), len(ptproduct)-1, array('d', ptproduct)),
            'ratio_n1': ROOT.TH1D(system + '_ratio_n1', ';p_{T,1}p_{T,2} (reco / gen);Probability', 80, .025, 2.025),
            'ratio_n2': ROOT.TH1D(system + '_ratio_n2', ';(p_{T,1}p_{T,2})^{2} (reco / gen);Probability', 80, .025, 2.025)}
        for hist in histograms[system].values():
            hist.Sumw2()
    for path in event_paths:
        with path.open() as stream:
            for line in stream:
                event = json.loads(line)
                if fold(event) != 1:
                    continue
                system = event['sample']
                for jet in event['jets']:
                    if not math.isclose(jet['radius'], .4, rel_tol=0., abs_tol=1e-6) or not 40 <= jet['pt'] < 80:
                        continue
                    if any('charge' not in track for population in ('reco', 'gen') for track in jet[population]):
                        raise ValueError('AN matching requires exported physical track charge')
                    stats[system]['jets'] += 1
                    for (a, ga), (b, gb) in matching_pairs(jet):
                        weight = event['generator_weight'] * event['event_weight'] * jet.get('jet_weight', 1.) * a['correction'] * b['correction']
                        reco, gen = a['pt'] * b['pt'], ga['pt'] * gb['pt']
                        histograms[system]['angular'].Fill(separation(a, b), separation(ga, gb), weight)
                        histograms[system]['momentum'].Fill(reco, gen, weight)
                        stats[system]['matched_pairs'] += 1
                        stats[system]['sum_weight'] += weight
                        for n in (1, 2):
                            ratio = (reco / gen) ** n
                            histograms[system]['ratio_n' + str(n)].Fill(ratio, weight)
                            stats[system]['sum_ratio_n' + str(n)] += weight * ratio
                            stats[system]['sum_ratio2_n' + str(n)] += weight * ratio * ratio
    directory = rootfile.mkdir('resolution')
    report = []
    for system, hist in histograms.items():
        if not stats[system]['jets']:
            continue
        directory.cd(); outdir = directory.mkdir(system); outdir.cd()
        drawer = ROOT.JDrawer(); ROOT.SetOwnership(drawer, False)
        drawer.SetLeftMargin(.18); drawer.SetRightMargin(.16); drawer.SetTitleOffsetY(1.3)
        drawer.SetHistogramAppearance(1., 1.3, .045, .045, .01, .001, .04, .04, 505, 505, 42)
        canvas = ROOT.TCanvas('resolution_' + system, '', 1200, 1000); canvas.Divide(2, 2)
        keep = []
        for panel, (name, original) in enumerate(hist.items(), 1):
            original.Write(name + '_weighted_counts')
            plot = original.Clone(name + '_normalized'); keep.append(plot)
            pad = canvas.cd(panel)
            pad.SetLeftMargin(.18); pad.SetRightMargin(.16); pad.SetBottomMargin(.18); pad.SetTopMargin(.11)
            if name in ('angular', 'momentum'):
                empty = 0
                for ix in range(1, plot.GetNbinsX() + 1):
                    total = sum(plot.GetBinContent(ix, iy) for iy in range(1, plot.GetNbinsY() + 1))
                    if total <= 0:
                        empty += 1
                        continue
                    for iy in range(1, plot.GetNbinsY() + 1):
                        plot.SetBinContent(ix, iy, plot.GetBinContent(ix, iy) / total)
                        plot.SetBinError(ix, iy, plot.GetBinError(ix, iy) / total)
                stats[system][name + '_unsupported_reco_rows'] = empty
                if name == 'angular':
                    plot.GetXaxis().SetRangeUser(.002, .4); plot.GetYaxis().SetRangeUser(.002, .4)
                else:
                    plot.GetXaxis().SetRangeUser(4, 400); plot.GetYaxis().SetRangeUser(4, 400)
                drawer.SetLogX(True); drawer.SetLogY(True)
                drawer.DrawHistogramToPad(plot, pad, original.GetXaxis().GetTitle(), original.GetYaxis().GetTitle(), system + ', 40-80 GeV; fixed-reco row normalized', 'COLZ')
            else:
                integral = plot.Integral()
                if integral > 0:
                    plot.Scale(1. / integral)
                drawer.SetLogX(False); drawer.SetLogY(False)
                drawer.DrawHistogramToPad(plot, pad, original.GetXaxis().GetTitle(), 'Probability', system + ', common reco jets; single-track corrected', 'HIST')
                n = name[-1]; total = stats[system]['sum_weight']
                mean = stats[system]['sum_ratio_n' + n] / total if total > 0 else None
                rms = math.sqrt(max(0., stats[system]['sum_ratio2_n' + n] / total - mean * mean)) if total > 0 else None
                stats[system]['mean_n' + n] = mean; stats[system]['rms_n' + n] = rms
                flow = original.GetBinContent(0) + original.GetBinContent(original.GetNbinsX() + 1)
                stats[system]['ratio_flow_weight_fraction_n' + n] = flow / total if total > 0 else None
                text = ROOT.TLatex(); text.SetNDC(); text.SetTextSize(.035)
                if mean is not None:
                    text.DrawLatex(.52, .82, f'mean = {mean:.4f}, RMS = {rms:.4f}')
                if total > 0 and flow > 0:
                    text.SetTextSize(.029)
                    text.DrawLatex(.23, .75, f'Outside plotted ratio range: {100*flow/total:.2f}% of weight')
                keep.append(text)
            heading = ROOT.TLatex(); heading.SetNDC(); heading.SetTextFont(42); heading.SetTextSize(.035)
            heading.DrawLatex(.18, .95, ('OO selected MB' if system == 'OO' else 'pp reference') + ', R=0.4, 40-80 GeV jets')
            heading.SetTextSize(.027)
            heading.DrawLatex(.18, .90, 'Single-track corrected; ' + ('fixed-reco rows normalized' if name in ('angular', 'momentum') else 'matched-pair momentum ratio'))
            keep.append(heading)
            plot.Write()
        filename = output / ('resolution_' + system + '_40_80.pdf')
        canvas.SaveAs(str(filename)); canvas.Write('canvas'); canvas.Close()
        report.append({'system': system, 'selection': 'selected MB' if system == 'OO' else 'ppRef', **stats[system], 'pdf': str(filename.resolve())})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='pair_counts.py JSON result')
    parser.add_argument('output', type=Path)
    parser.add_argument('--jet-bin', type=int, choices=(0, 1, 2), help='first 40-80 closure: use 0; omit for all')
    parser.add_argument('--events', type=Path, nargs='*', default=[], help='native exported JSONL for compact AN matching diagnostics')
    parser.add_argument('--tables-only', action='store_true')
    args = parser.parse_args()
    result = json.loads(args.input.read_text())
    # A per-system streaming derivation retains the standard empty grid for API
    # compatibility; do not publish plots for a system absent from its inputs.
    systems = {row['key'][0] for row in result['map'] + result['jet_independent_map']}
    result['holdout_spectra'] = [row for row in result['holdout_spectra'] if row['key'][0] in systems]
    args.output.mkdir(parents=True, exist_ok=True)
    ROOT = root_setup()
    an_source = Path('/eos/user/j/jvelazqu/HIN/Plan/AN2023_081_v13.pdf')
    report = {'input': str(args.input.resolve()), 'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
              'plotter_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'uncertainties': 'physical-event statistical bars conditional on fixed map and common jet normalization; finite-map uncertainty not evaluated',
              'open_marker_policy': 'open squares are supported pair contributions only, not full corrected EEC closure',
              'AN_reference': 'AN2023_081_v13.pdf Sec.7.5.1 (printed pp.56-65), source /eos/user/j/jvelazqu/HIN/Plan/AN2023_081_v13.pdf; matching steps1-7 lines788-794; existing src/EECAnalyzer.cxx ConstructParticleResponses and plotting/trackParticleMatchingPlotter.C',
              'AN_sha256': hashlib.sha256(an_source.read_bytes()).hexdigest() if an_source.exists() else None,
              'shape_convention': 'AN Sec7.5.1 line810 and Sec7.1 line439; compareEECinDefinedBins.C whole-bin normalization selected by drawingRange(.008,.39), effective stored-edge interval [.008,.4)',
              'tables': write_tables(ROOT, result, args.output)}
    if not args.tables_only:
        rootfile = ROOT.TFile(str(args.output / 'closure.root'), 'RECREATE')
        report['closure'] = draw_closure(ROOT, result, args.output, rootfile, args.jet_bin)
        shape_result, shape_statistics, shape_support = normalized_shapes(result)
        report['shape_closure'] = draw_closure(ROOT, shape_result, args.output, rootfile, args.jet_bin, shape=True, statistics_override=shape_statistics)
        report['shape_support'] = shape_support
        report['shape_uncertainties'] = 'Joint physical-event delta-method covariance includes each curve normalization and numerator/truth covariance; conditional on fixed training map'
        if args.events:
            report['resolution'] = resolution(ROOT, result, args.events, args.output, rootfile)
            report['resolution_method'] = 'Existing greedy matching: descending reco pT; unmatched equal-charge gen within dR<0.05 and |dpT|<0.5*pT_reco; select minimum |dpT|. Independent validation events, common selected reco jets 40-80 GeV, single-track factor once. Empty response rows unfilled; ratios n=1,2; no new fit.'
        rootfile.cd(); put_json(ROOT, 'Provenance', report); rootfile.Close()
    (args.output / 'support_and_closure.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'tables': len(report['tables']), 'closure_plots': len(report.get('closure', [])),
                      'shape_plots': len(report.get('shape_closure', [])), 'resolution_plots': len(report.get('resolution', [])), 'output': str(args.output.resolve())}))


if __name__ == '__main__':
    main()
