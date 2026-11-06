#!/usr/bin/env python3
"""Project the standard EEC manager and plot bounded OO GenGen closure.

Run from the repository root, in its configured ROOT environment. No physics
analysis is reimplemented: central values use EECHistogramManager projections;
the opt-in analyzer ledgers supply joint physical-event statistical moments.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from array import array
from collections import Counter, defaultdict
from pathlib import Path

from closure_statistics import (ALL, COEFFICIENTS, JET_BINS, METRICS, SELECTIONS,
                                UNCERTAINTY, cluster_vectors, comparison, covariance, derive,
                                joint_moments)

COMPARISONS = (("Bmix", "Btrue"), ("Bperp", "Btrue"), ("Smix", "PP"),
               ("Sperp", "PP"), ("mix_PMavg", "PH"),
               ("mix_same_avg", "HH"), ("perp_PMavg", "PH"),
               ("perp_same_avg", "HH"), ("PP", "T"), ("PH", "T"), ("HH", "T"))
PAIRS = (0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
SELECTION_LABELS = dict(selected_integrated='selected-event integrated',
                        C0_30='0-30% centrality', C30_50='30-50% centrality',
                        C50_90='50-90% centrality')


def rows(path):
    with path.open(newline="") as stream:
        yield from csv.DictReader(stream)


def write_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def jet_bin(value, bins):
    return next((i for i, (low, high) in enumerate(bins) if low <= value < high), None)


def selection_names(centrality, selections):
    return [name for name, low, high in selections if low <= centrality < high]


def read_ledgers(prefix, bins, selections, first_track_bin):
    partners = list(rows(Path(str(prefix) + ".partners.csv")))
    pool = list(rows(Path(str(prefix) + ".pool.csv")))
    jets = list(rows(Path(str(prefix) + ".jets.csv")))
    partner_keys = [row["signal_key"] for row in partners]
    if len(set(partner_keys)) != len(partner_keys):
        raise ValueError("Duplicate physical signal identity in partner ledger")
    # Audit the actual saved assignments, preserving decimal identities as
    # strings/integers. Source paths and file entries never establish independence.
    signal_physical = {(row['run'], row['lumi'], row['event']) for row in partners
                       if all(field in row for field in ('run', 'lumi', 'event'))}
    known_signal_backgrounds = {row['signal_background_key'] for row in partners
        if int(row.get('signal_background_signature', '0:0:0').split(':', 1)[0]) > 0}
    for row in partners:
        if row.get('status') != 'matched' or 'donor1_run' not in row:
            continue
        donor_ids = [tuple(row['donor{}_{}'.format(index, field)] for field in ('run', 'lumi', 'event'))
                     for index in (1, 2)]
        if donor_ids[0] == donor_ids[1] or any(identity in signal_physical for identity in donor_ids):
            raise ValueError('Saved partner assignment contains self/duplicate/signal-overlap physical identities')
        donor_backgrounds = [row.get('donor{}_background_key'.format(index), '0') for index in (1, 2)]
        known = [key for key in donor_backgrounds if key != '0']
        if len(known) != len(set(known)) or any(key in known_signal_backgrounds for key in known):
            raise ValueError('Saved partner assignment contains observed shared signal/donor background')
    cluster_keys, background_keys = {}, {}
    for row in partners:
        identity = row['signal_key']
        background = row.get('signal_background_key', '0')
        count = int(row.get('signal_background_signature', '0:0:0').split(':', 1)[0])
        background_keys[identity] = background
        cluster_keys[identity] = ('background:' + background if count > 0
                                  else 'physical:' + identity)
    jet_lookup, exposure = {}, defaultdict(lambda: dict(selected_jets=0, matched_jets=0,
                                                        selected_events=set(), matched_events=set(),
                                                        selected_sumw=0., matched_sumw=0.))
    for row in jets:
        key = row["signal_key"], int(row["jet_index"])
        if key in jet_lookup:
            raise ValueError("Duplicate signal physical-event / jet identity")
        jet_lookup[key] = row
        ibin = jet_bin(float(row["jet_pt"]), bins)
        if ibin is None:
            continue
        for selection in selection_names(float(row["centrality"]), selections):
            current = exposure[selection, ibin]
            current["selected_jets"] += 1
            current["selected_events"].add(key[0])
            current["selected_sumw"] += float(row["event_weight"])
            if row["status"] == "matched":
                current["matched_jets"] += 1
                current["matched_events"].add(key[0])
                current["matched_sumw"] += float(row["event_weight"])
    raw_events, event_cent = defaultdict(lambda: defaultdict(float)), {}
    source_rows = 0
    for row in rows(Path(str(prefix) + ".contributions.csv")):
        source_rows += 1
        identity = row["signal_key"]
        ledger_jet = jet_lookup.get((identity, int(row["jet_index"])))
        if ledger_jet is None or ledger_jet["status"] != "matched":
            raise ValueError("Contribution outside common successfully matched jet sample")
        if int(row.get("track_pt_bin", first_track_bin)) < first_track_bin:
            continue
        ibin = jet_bin(float(row["jet_pt"]), bins)
        if ibin is None:
            continue
        centrality = float(row["centrality"])
        if identity in event_cent and event_cent[identity] != centrality:
            raise ValueError("One physical event has inconsistent centrality")
        event_cent[identity] = centrality
        origin = int(row["subevent"])
        if origin not in range(4):
            raise ValueError("Ledger must retain raw ordered subevent categories 0..3")
        key = identity, ibin, int(row["power"]), int(row["delta_r_bin"])
        raw_events[key][int(row["pairing"]), origin] += float(row["value"])
    event_vectors = defaultdict(dict)
    for (identity, ibin, power, drbin), raw in raw_events.items():
        values = derive(raw)
        for selection in selection_names(event_cent[identity], selections):
            event_vectors[selection, ibin, power, drbin][identity] = values
    for current in exposure.values():
        for key in ("selected_events", "matched_events"):
            current[key] = len(current[key])
        current["matching_coverage"] = (current["matched_jets"] / current["selected_jets"]
                                        if current["selected_jets"] else None)
    return dict(partners=partners, pool=pool, jets=jets, exposure=exposure,
                event_vectors=event_vectors, contribution_rows=source_rows,
                cluster_keys=cluster_keys, background_keys=background_keys)


def save_postprocessing_snapshot(output, ROOT):
    """Freeze the actual postprocessor separately from the native-run snapshot."""
    import matplotlib
    root = Path(__file__).resolve().parents[2]
    snapshot = output / 'postprocessing_sources'
    sources = ('validation/run3Mixing/plot_closure_pilot.py',
               'validation/run3Mixing/closure_statistics.py',
               'validation/run3Mixing/test_closure_statistics.py',
               'validation/run3Mixing/test_closure_ledgers.py',
               'validation/run3Mixing/test_raw_manager_projection.py',
               'plotting/EECHistogramManager.h', 'plotting/EECHistogramManager.cxx')
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                    source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
                    working_directory=str(Path.cwd()), command=shlex.join([sys.executable] + sys.argv),
                    argv=[sys.executable] + sys.argv, python=platform.python_version(),
                    ROOT=ROOT.gROOT.GetVersion(), matplotlib=matplotlib.__version__, files=[])
    for relative in sources:
        source, target = root / relative, snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        manifest['files'].append(dict(path=relative, snapshot=str(target.relative_to(output)),
                                      sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    library = root / 'plotting/DrawingClasses.so'
    manifest['loaded_plotting_library'] = dict(path=str(library), sha256=hashlib.sha256(library.read_bytes()).hexdigest())
    write_json(output / 'postprocessing_manifest.json', manifest)
    with (output / 'commands.sh').open('x') as stream:
        stream.write('#!/usr/bin/env bash\nset -euo pipefail\n')
        stream.write('# Actual executed command; choose a fresh --output-dir when reproducing.\n')
        stream.write('cd ' + shlex.quote(str(Path.cwd())) + '\n' + manifest['command'] + '\n')
    return manifest


def initialize_root():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if ROOT.gSystem.Load("plotting/DrawingClasses.so") < 0:
        raise RuntimeError("Build plotting/DrawingClasses.so before plotting")
    ROOT.gInterpreter.Declare('#include "plotting/EECHistogramManager.h"\n#include "plotting/EECCard.h"')
    return ROOT


def project(ROOT, source, bins, selections, powers, threshold):
    """Use manager projections, summing disjoint native centrality intervals."""
    projected, retained = {}, []
    card = ROOT.EECCard(source)
    if "OO" not in str(card.GetDataType()) or "MC" not in str(card.GetDataType()):
        raise ValueError("Requires OO MC")
    correlation = card.GetCardEntry(ROOT.EECCard.kMcCorrelationType)
    if not correlation or round(correlation.GetMatrixArray()[0]) != 3:
        raise ValueError("Requires native GenGen jets and generator particles")
    if card.GetBackgroundMethods() & 6 != 6:
        raise ValueError("Requires mixed and perpendicular categories in the same ROOT file")
    if abs(card.GetJetRadius() - .4) > 1e-5:
        raise ValueError("This pilot requires R=0.4")
    track = card.GetBinIndexTrackPtEEC(threshold)
    if track < 0 or abs(card.GetLowBinBorderTrackPtEEC(track) - threshold) > 1e-5:
        raise ValueError("Requested exact pair-pT threshold unavailable")
    native = [(card.GetLowBinBorderCentrality(i), card.GetHighBinBorderCentrality(i))
              for i in range(card.GetNCentralityBins())]
    jet_indices = [card.FindBinIndexJetPtEEC((low, high)) for low, high in bins]
    if any(index < 0 for index in jet_indices):
        raise ValueError("Requested jet bins must exist exactly in the stored card")
    for power in powers:
        if card.FindWeightExponentIndex(power) < 1:
            raise ValueError("Requested power unavailable")
        current_card = ROOT.EECCard(source)
        manager = ROOT.EECHistogramManager(source, current_card)
        manager.SetLoadEnergyEnergyCorrelators(True)
        manager.SetLoadedWeightExponent(power)
        # Set indices from the native stored axes through the standard API.
        for setter in (manager.SetCentralityBins, manager.SetJetPtBinsEEC,
                       manager.SetTrackPtBinsEEC):
            setter(True, 0, array('d', [0.]), True)
        manager.SetCentralityBinRange(0, current_card.GetNCentralityBins() - 1)
        manager.SetJetPtBinRangeEEC(min(jet_indices), max(jet_indices))
        manager.SetTrackPtBinRangeEEC(track, track)
        manager.LoadRawEnergyEnergyCorrelators()
        retained.extend((current_card, manager))
        for selection, low, high in selections:
            selected = [i for i, (left, right) in enumerate(native)
                        if left >= low - 1e-6 and right <= high + 1e-6]
            if not selected or abs(native[selected[0]][0] - low) > 1e-6 or abs(native[selected[-1]][1] - high) > 1e-6:
                raise ValueError("Selection {} does not align with stored centrality intervals".format(selection))
            for ibin, jet_index in enumerate(jet_indices):
                categories = {}
                for pairing in PAIRS:
                    for origin in (0, 1, 3, ALL):
                        result = None
                        for cent in selected:
                            histogram = manager.GetHistogramEnergyEnergyCorrelator(0, cent, jet_index,
                                                                                 track, pairing, origin)
                            if not histogram:
                                raise ValueError("Missing pairing {} origin {}".format(pairing, origin))
                            if result is None:
                                result = histogram.Clone("{}_J{}_n{}_pair{}_S{}".format(selection, ibin, power, pairing, origin))
                                result.SetDirectory(0)
                            else:
                                result.Add(histogram)
                        categories[pairing, origin] = result
                metrics = {}
                for name, mapping in COEFFICIENTS.items():
                    histogram = categories[0, ALL].Clone("{}_J{}_n{}_{}".format(selection, ibin, power, name))
                    histogram.SetDirectory(0)
                    histogram.Reset("ICES")
                    for key, coefficient in mapping.items():
                        histogram.Add(categories[key], coefficient)
                    metrics[name] = histogram
                projected[selection, ibin, power] = categories, metrics
    return projected, retained


def build_summary(projected, ledgers, bins, powers, selections, output, ROOT, source_name, provenance=None):
    summary = dict(status="BOUNDED_MC_DIAGNOSTIC_NOT_FROZEN", input_root=str(source_name),
                   postprocessing=provenance,
                   uncertainty_method=UNCERTAINTY, normalization="raw shared signal-jet exposure; no independent component normalization",
                   origin_convention="manager S1 already merges raw PH and HP; PH counted once",
                   covariance_scope="Joint metric covariance per Delta-r bin and integrated interval; full cross-bin/power/jet matrix NOT_RUN; retained event vectors and cluster keys permit reconstruction",
                   systematics_join="joint_event_vectors.jsonl signal_key joins original .partners.csv donor identities and reuse; covariance_cluster_key joins observed signal background clusters",
                   pair_threshold_GeV=2., centrality_note="Selected-event acceptance only; no unsupported 0-100% physical coverage claim",
                   empty_bins="Retained as zero yields; ratios to nonpositive truth are null and explicitly masked",
                   ledger_contribution_rows=ledgers["contribution_rows"],
                   raw_pair_sumw2_errors="Diagnostic only; replaced by joint physical-signal-event moments in plotted derived curves",
                   channels=[])
    root_output = ROOT.TFile.Open(str(output / "closure.root"), "CREATE")
    if not root_output or root_output.IsZombie():
        raise RuntimeError("Cannot create derived ROOT output")
    ROOT.TNamed("uncertainty_method", UNCERTAINTY).Write()
    per_bin_rows, integral_rows = [], []
    with (output / "joint_event_vectors.jsonl").open("x") as stream:
        for key, events in sorted(ledgers["event_vectors"].items()):
            for identity, values in sorted(events.items()):
                stream.write(json.dumps(dict(signal_key=identity,
                                             signal_background_key=ledgers['background_keys'].get(identity, '0'),
                                             covariance_cluster_key=ledgers['cluster_keys'].get(identity, 'physical:' + identity),
                                             selection=key[0], jet_bin=key[1], power=key[2],
                                             delta_r_bin=key[3], values=values), allow_nan=False) + "\n")
    for (selection, ibin, power), (categories, metrics) in projected.items():
        template = metrics["T"]
        nbin = template.GetNbinsX()
        current = dict(selection=selection, jet_pt_low=bins[ibin][0], jet_pt_high=bins[ibin][1], power=power,
                       exposure=ledgers["exposure"].get((selection, ibin), dict(selected_jets=0, matched_jets=0,
                          selected_events=0, matched_events=0, selected_sumw=0., matched_sumw=0., matching_coverage=None)),
                       bins=[], integrals={}, maximum_manager_ledger_absolute_difference=0.)
        scope = "{}_J{}-{}_n{}".format(selection, int(bins[ibin][0]), int(bins[ibin][1]), power)
        directory = root_output.mkdir(scope)
        directory.cd()
        exposure_histogram = ROOT.TH1D(scope + '_exposure', 'Common signal-jet exposure', 6, .5, 6.5)
        for index, field in enumerate(('selected_jets', 'matched_jets', 'selected_events', 'matched_events',
                                       'selected_sumw', 'matched_sumw'), 1):
            exposure_histogram.GetXaxis().SetBinLabel(index, field)
            exposure_histogram.SetBinContent(index, current['exposure'][field])
        exposure_histogram.Write()
        for histogram in categories.values():
            histogram.Write()
        masks = {}
        for estimate, truth in COMPARISONS:
            mask = template.Clone("{}_{}_over_{}_invalid_mask".format(scope, estimate, truth))
            mask.Reset("ICES"); mask.SetDirectory(0)
            masks[estimate, truth] = mask
        for drbin in range(nbin + 2):
            events = ledgers["event_vectors"].get((selection, ibin, power, drbin), {})
            clustered = cluster_vectors(events, ledgers['cluster_keys'])
            total, cov = joint_moments(clustered.values())
            width = template.GetBinWidth(drbin)
            # Independent ledger/manager agreement catches subevent double counting,
            # threshold mistakes, and accidental exposure differences.
            for name, histogram in metrics.items():
                stored = histogram.GetBinContent(drbin) * width
                scale = sum(abs(coefficient * categories[key].GetBinContent(drbin) * width)
                            for key, coefficient in COEFFICIENTS[name].items())
                current['maximum_manager_ledger_absolute_difference'] = max(
                    current['maximum_manager_ledger_absolute_difference'], abs(stored - total[name]))
                if abs(stored - total[name]) > max(1e-5, 3e-5 * scale):
                    raise ValueError("Manager/ledger disagreement {} {} bin{}: {} versus {}".format(scope, name, drbin, stored, total[name]))
                histogram.SetBinError(drbin, math.sqrt(max(0., covariance(cov, name, name))) / width)
            total = {name: histogram.GetBinContent(drbin) * width for name, histogram in metrics.items()}
            comparisons = [comparison(total, cov, estimate, truth) for estimate, truth in COMPARISONS]
            for check in comparisons:
                masks[check["estimate"], check["truth"]].SetBinContent(drbin, 0 if check["valid_ratio"] else 1)
            entry = dict(bin=drbin, low=template.GetXaxis().GetBinLowEdge(drbin),
                         high=template.GetXaxis().GetBinUpEdge(drbin), underflow=drbin == 0, overflow=drbin == nbin + 1,
                         raw_yields=total, raw_yield_errors={name: math.sqrt(max(0., covariance(cov, name, name))) for name in METRICS},
                         comparisons=comparisons,
                         covariance_upper=[[first, second, value] for (first, second), value in cov.items() if value])
            current["bins"].append(entry)
            row = dict(selection=selection, jet_pt_low=bins[ibin][0], jet_pt_high=bins[ibin][1], power=power,
                       delta_r_bin=drbin, delta_r_low=entry['low'], delta_r_high=entry['high'])
            for name in METRICS:
                row[name] = total[name]
                row[name + "_error"] = entry["raw_yield_errors"][name]
            for check in comparisons:
                for name in ("ratio", "ratio_error", "residual", "residual_error", "valid_ratio"):
                    row[check['estimate'] + "_" + name] = check[name]
            per_bin_rows.append(row)
        for view, low, high in (("full", template.GetXaxis().GetXmin(), template.GetXaxis().GetXmax()),
                                ("legacy", .008, .39)):
            first = max(1, template.FindFixBin(low + 1e-10))
            last = min(nbin, template.FindFixBin(high - 1e-10))
            integrated = defaultdict(lambda: defaultdict(float))
            for drbin in range(first, last + 1):
                for identity, vector in ledgers["event_vectors"].get((selection, ibin, power, drbin), {}).items():
                    for name, value in vector.items():
                        integrated[identity][name] += value
            clustered = cluster_vectors(integrated, ledgers['cluster_keys'])
            total, cov = joint_moments(clustered.values())
            total = {name: histogram.Integral(first, last, 'width') for name, histogram in metrics.items()}
            checks = [comparison(total, cov, estimate, truth) for estimate, truth in COMPARISONS]
            current['integrals'][view] = dict(requested_delta_r=[low, high],
                actual_bin_interval=[template.GetXaxis().GetBinLowEdge(first), template.GetXaxis().GetBinUpEdge(last)],
                bin_indices=[first, last], raw_yields=total,
                raw_yield_errors={name: math.sqrt(max(0., covariance(cov, name, name))) for name in METRICS},
                comparisons=checks, covariance_upper=[[a, b, v] for (a, b), v in cov.items() if v])
            for name in METRICS:
                row = dict(selection=selection, jet_pt_low=bins[ibin][0], jet_pt_high=bins[ibin][1], power=power,
                           view=view, metric=name, integral=total[name],
                           integral_error=math.sqrt(max(0., covariance(cov, name, name))),
                           matched_jets=current['exposure']['matched_jets'], selected_jets=current['exposure']['selected_jets'],
                           matching_coverage=current['exposure']['matching_coverage'], uncertainty_method=UNCERTAINTY)
                match = next((check for check in checks if check['estimate'] == name), None)
                for field in ('truth', 'ratio', 'ratio_error', 'residual', 'residual_error', 'valid_ratio'):
                    row[field] = match[field] if match else None
                integral_rows.append(row)
        for histogram in list(metrics.values()) + list(masks.values()):
            histogram.Write()
        current['status'] = "MATCHED_JETS" if current['exposure']['matched_jets'] else "EMPTY_NO_MATCHED_JETS"
        regular = [entry for entry in current['bins'] if not entry['underflow'] and not entry['overflow']]
        current['positive_background_truth_bins'] = sum(entry['raw_yields']['Btrue'] > 0 for entry in regular)
        current['background_status'] = ('SUPPORTED_DIAGNOSTIC_ONLY' if current['positive_background_truth_bins']
                                        else 'EMPTY_TRUTH_BACKGROUND_NO_CLOSURE_TEST')
        current['invalid_ratio_bin_counts'] = {estimate: sum(not next(check for check in entry['comparisons']
                if check['estimate'] == estimate)['valid_ratio'] for entry in regular)
                for estimate, _ in COMPARISONS}
        summary['channels'].append(current)
    root_output.Close()
    for filename, contents in (("bin_summary.csv", per_bin_rows), ("closure_summary.csv", integral_rows)):
        with (output / filename).open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(contents[0]))
            writer.writeheader(); writer.writerows(contents)
    write_json(output / 'summary.json', summary)
    return summary


def draw_groups(summary, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = ('black', '#2166ac', '#b2182b', '#1b7837')
    for channel in summary['channels']:
        regular = [entry for entry in channel['bins'] if not entry['underflow'] and not entry['overflow']]
        for view in ('legacy', 'full'):
            low, high = channel['integrals'][view]['requested_delta_r']
            visible = [entry for entry in regular if entry['high'] > low and entry['low'] < high]
            x = [math.sqrt(max(entry['low'], 1e-12) * entry['high']) if entry['low'] > 0 else .5 * entry['high'] for entry in visible]
            width = [entry['high'] - entry['low'] for entry in visible]
            groups = [('background', ('Btrue', 'Bmix', 'Bperp'), ('Truth PH+HH', 'Mixed cones', 'Perpendicular cones')),
                      ('signal', ('PP', 'Smix', 'Sperp'), ('Truth PP', 'T - mixed background', 'T - perpendicular background')),
                      ('decomposition', ('T', 'PP', 'PH', 'HH'), ('Total T', 'Pythia-Pythia (PP)', 'Pythia-Hijing (PH, counted once)', 'Hijing-Hijing (HH)'))]
            for group, names, labels in groups:
                figure, (upper, lower) = plt.subplots(2, 1, figsize=(7.6, 7.4), sharex=True,
                    gridspec_kw={'height_ratios': [3., 1.25]}, constrained_layout=True)
                for color, name, label in zip(colors, names, labels):
                    y = [entry['raw_yields'][name] / w for entry, w in zip(visible, width)]
                    errors = [entry['raw_yield_errors'][name] / w for entry, w in zip(visible, width)]
                    upper.errorbar(x, y, yerr=errors, color=color, marker='o', markersize=2.5,
                                   linewidth=.8, elinewidth=.6, label=label)
                upper.axhline(0., color='0.6', linewidth=.6)
                upper.set_ylabel('Weighted pair yield / $\\Delta r$\n(common matched signal-jet exposure)')
                upper.legend(fontsize=8)
                upper.set_title('OO GenGen R=0.4; {}\n{} < jet $p_T$ < {} GeV; $p_T^{{ch}}>2$ GeV; n={}\n{} matched / {} selected jets; conditional fixed-pool errors'.format(
                    SELECTION_LABELS.get(channel['selection'], channel['selection']), int(channel['jet_pt_low']), int(channel['jet_pt_high']), channel['power'],
                    channel['exposure']['matched_jets'], channel['exposure']['selected_jets']), fontsize=10)
                for color, name, label in zip(colors[1:], names[1:], labels[1:]):
                    checks = [next(check for check in entry['comparisons'] if check['estimate'] == name) for entry in visible]
                    valid = [(xx, check) for xx, check in zip(x, checks) if check['valid_ratio']]
                    if valid:
                        lower.errorbar([item[0] for item in valid], [item[1]['ratio'] for item in valid],
                            yerr=[item[1]['ratio_error'] for item in valid], color=color, marker='o',
                            markersize=2.5, linewidth=.8, elinewidth=.6)
                for entry in visible:
                    truth = entry['raw_yields'][names[0]]
                    if truth <= 0:
                        lower.axvspan(max(entry['low'], max(low, 1e-6)), entry['high'], color='0.8', alpha=.5)
                lower.axhline(1., color='0.4', linestyle='--', linewidth=.8)
                lower.set_ylabel('Fraction of T' if group == 'decomposition' else 'Estimate / truth')
                lower.set_xlabel('$\\Delta r$; grey = invalid ratio (truth $\\leq 0$)')
                if not any(entry['raw_yields'][names[0]] > 0 for entry in visible):
                    lower.text(.5, .5, 'No valid truth denominator', transform=lower.transAxes, ha='center')
                if not channel['exposure']['matched_jets']:
                    upper.text(.5, .5, 'EMPTY: no successfully matched jets', transform=upper.transAxes, ha='center')
                if low > 0:
                    lower.set_xscale('log')
                lower.set_xlim(low if low > 0 else 0., high)
                lower.grid(axis='y', alpha=.25)
                stem = '{}_J{}-{}_n{}_{}_{}'.format(channel['selection'], int(channel['jet_pt_low']), int(channel['jet_pt_high']), channel['power'], group, view)
                for extension in ('pdf', 'png'):
                    figure.savefig(output / (stem + '.' + extension), dpi=150)
                plt.close(figure)
            figure, axes = plt.subplots(3, 1, figsize=(7.6, 9.), sharex=True, constrained_layout=True)
            component_groups = ((('PH', 'mix_PMavg', 'perp_PMavg'), ('Truth PH', 'Average PM1,PM2', 'Average PP1,PP2')),
                                (('mix_HMavg_minus_cross', 'perp_HMavg_minus_cross'), ('Average HM - M1M2', 'Average HP - P1P2')),
                                (('HH', 'mix_same_avg', 'perp_same_avg'), ('Truth HH', 'Average M1M1,M2M2', 'Average P1P1,P2P2')))
            for axis, (names, labels) in zip(axes, component_groups):
                current_colors = colors if len(names) == 3 else colors[1:]
                for color, name, label in zip(current_colors, names, labels):
                    axis.errorbar(x, [entry['raw_yields'][name] / w for entry, w in zip(visible, width)],
                        yerr=[entry['raw_yield_errors'][name] / w for entry, w in zip(visible, width)],
                        color=color, marker='o', markersize=2.5, linewidth=.8, elinewidth=.6, label=label)
                axis.axhline(0., color='0.5', linewidth=.6)
                axis.legend(fontsize=8); axis.set_ylabel('Weighted pair yield / $\\Delta r$')
            axes[0].set_title('OO GenGen R=0.4; {}; J{}-{} GeV; n={}\nCommon matched exposure; conditional fixed-pool errors; signed values retained'.format(
                SELECTION_LABELS.get(channel['selection'], channel['selection']), int(channel['jet_pt_low']), int(channel['jet_pt_high']), channel['power']), fontsize=10)
            axes[-1].set_xlim(low, high)
            if low > 0:
                axes[-1].set_xscale('log')
            axes[-1].set_xlabel('$\\Delta r$')
            stem = '{}_J{}-{}_n{}_components_{}'.format(channel['selection'], int(channel['jet_pt_low']), int(channel['jet_pt_high']), channel['power'], view)
            for extension in ('pdf', 'png'):
                figure.savefig(output / (stem + '.' + extension), dpi=150)
            plt.close(figure)


def draw_qa(ledgers, summary, output, bins):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    partners, pool, jets = (ledgers[key] for key in ('partners', 'pool', 'jets'))
    matched = [row for row in partners if row['status'] == 'matched']
    reasons = Counter(row['reason'] for row in partners if row['status'] != 'matched')
    qa = dict(signal_events=len(partners), matched_signal_events=len(matched), failures=dict(reasons),
              observed_signal_background_clusters=len(set(ledgers['cluster_keys'].values())),
              observed_matched_signal_background_clusters=len({ledgers['cluster_keys'][row['signal_key']] for row in matched}),
              pool_entries=len(pool), unique_used_donors=sum(int(row['uses']) > 0 for row in pool),
              donor_uses=sum(int(row['uses']) for row in pool), max_reuse=max([int(row['uses']) for row in pool] or [0]),
              donor_pool_weight_note="Donor generator weights retained in ledger; not multiplied into signal-jet pair terms",
              shared_embedding_background_provenance="UNKNOWN: no verified shared-underlying-HIJING mapping; independence not certified",
              event_identity_note="Physical run/lumi/event retained as integer/string identities; file entry is locator only",
              uncertainty_method=UNCERTAINTY)
    qa['signal_centrality_acceptance'] = {}
    for status, label in ((None, 'selected'), ('matched', 'matched')):
        sample = [row for row in partners if status is None or row['status'] == status]
        qa['signal_centrality_acceptance'][label] = dict(
            hiBin_0_to_179=sum(0 <= int(row['hiBin']) < 180 for row in sample),
            hiBin_180_to_199=sum(180 <= int(row['hiBin']) < 200 for row in sample),
            invalid_hiBin=sum(not 0 <= int(row['hiBin']) < 200 for row in sample))
    qa['pool_occupancy'] = {}
    for fold in ('A', 'B'):
        sample = [row for row in pool if row['pool'] == fold]
        qa['pool_occupancy'][fold] = dict(entries=len(sample),
            eligible=sum(row['eligible'] == '1' for row in sample),
            unique_used=sum(int(row['uses']) > 0 for row in sample), uses=sum(int(row['uses']) for row in sample))
    qa['pool_generator_weights'] = dict(
        minimum=min([float(row['generator_weight']) for row in pool] or [0.]),
        maximum=max([float(row['generator_weight']) for row in pool] or [0.]),
        nonunity=sum(float(row['generator_weight']) != 1. for row in pool))
    nh_values = [float(row['nh']) for row in partners]
    nh_low, nh_high = min(nh_values or [0.]), max(nh_values or [1.])
    nh_width = max(1., (nh_high - nh_low + 1.) / 20.)
    activity = []
    for activity_bin in range(20):
        left, right = nh_low + activity_bin * nh_width, nh_low + (activity_bin + 1) * nh_width
        signals = [row for row in partners if left <= float(row['nh']) < right]
        matched_count = sum(row['status'] == 'matched' for row in signals)
        channel_rows = []
        signal_keys = {row['signal_key'] for row in signals}
        for ibin, (low, high) in enumerate(bins):
            selected = [row for row in jets if row['signal_key'] in signal_keys and low <= float(row['jet_pt']) < high]
            passed = sum(row['status'] == 'matched' for row in selected)
            channel_rows.append(dict(jet_pt=[low, high], selected_jets=len(selected), matched_jets=passed,
                                     coverage=passed / len(selected) if selected else None))
        activity.append(dict(nh_low=left, nh_high=right, selected_signals=len(signals), matched_signals=matched_count,
                             coverage=matched_count / len(signals) if signals else None, jet_bins=channel_rows))
    qa['matching_acceptance_vs_activity'] = activity
    write_json(output / 'matching_qa.json', qa)
    figure, axes = plt.subplots(2, 3, figsize=(14., 8.), constrained_layout=True)
    axes = axes.ravel()
    for partner in (1, 2):
        axes[0].hist([float(row['donor{}_vz'.format(partner)]) - float(row['vz']) for row in matched],
                     bins=30, histtype='step', label='Donor {}'.format(partner))
        axes[1].hist([float(row['donor{}_nh'.format(partner)]) - float(row['nh']) for row in matched],
                     bins=30, histtype='step', label='Donor {}'.format(partner))
    axes[0].set_xlabel('$v_z^{donor}-v_z^{signal}$ [cm]'); axes[0].legend()
    axes[1].set_xlabel('$N_H^{donor}-N_H^{signal}$'); axes[1].legend()
    axes[2].hist([int(row['uses']) for row in pool], bins=range(0, qa['max_reuse'] + 2), histtype='step')
    axes[2].set_xlabel('Uses per physical donor (zero included)')
    for status, label in ((None, 'Selected signals'), ('matched', 'Matched signals')):
        axes[3].hist([float(row['nh']) for row in partners if status is None or row['status'] == status],
                     bins=30, histtype='step', label=label)
    axes[3].set_xlabel('Stable charged HIJING $N_H$; $p_T>0.7$, $|\\eta|<2.4$'); axes[3].legend()
    integrated = next(selection for selection, _, _ in SELECTIONS if selection == 'selected_integrated')
    labels, coverage = [], []
    for ibin, (low, high) in enumerate(bins):
        current = ledgers['exposure'].get((integrated, ibin), {})
        labels.append('{}-{}'.format(int(low), int(high)))
        coverage.append(current.get('matching_coverage') or 0.)
    axes[4].bar(labels, coverage); axes[4].set_ylim(0., 1.)
    axes[4].set_xlabel('Jet $p_T$ [GeV]'); axes[4].set_ylabel('Matched / selected jets')
    for status, label in ((None, 'Selected signals'), ('matched', 'Matched signals')):
        axes[5].hist([float(row['hiBin']) / 2. for row in partners if status is None or row['status'] == status],
                     bins=[-1., 0., 30., 50., 90., 100.], histtype='step', label=label)
    axes[5].set_xlabel('Stored signal centrality (hiBin/2)'); axes[5].legend()
    figure.suptitle('Bounded OO GenGen mixing QA; {} / {} matched signals; {} used physical donors\nUnknown shared embedding-background provenance; finite matching tolerances'.format(
        len(matched), len(partners), qa['unique_used_donors']))
    for extension in ('pdf', 'png'):
        figure.savefig(output / ('matching_QA.' + extension), dpi=150)
    plt.close(figure)
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--ledger-prefix', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--jet-bins', default='40:80,80:120,120:160')
    parser.add_argument('--powers', default='1,2')
    parser.add_argument('--pair-threshold', type=float, default=2.)
    args = parser.parse_args()
    bins = tuple(tuple(float(value) for value in interval.split(':')) for interval in args.jet_bins.split(','))
    powers = tuple(int(value) for value in args.powers.split(','))
    if args.pair_threshold != 2.:
        raise ValueError('Pilot ledger/plot contract currently requires pTch>2 GeV')
    ROOT = initialize_root()
    source = ROOT.TFile.Open(str(args.input), 'READ')
    if not source or source.IsZombie():
        raise ValueError('Cannot read input ROOT file: ' + str(args.input))
    sparse = source.Get('energyEnergyCorrelator')
    if not sparse:
        raise ValueError('Input must be raw analyzer output with standard EEC THnSparse')
    first_track = sparse.GetAxis(2).FindFixBin(args.pair_threshold + 1e-6)
    ledgers = read_ledgers(args.ledger_prefix, bins, SELECTIONS, first_track)
    projected, retained = project(ROOT, source, bins, SELECTIONS, powers, args.pair_threshold)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    provenance = save_postprocessing_snapshot(args.output_dir, ROOT)
    summary = build_summary(projected, ledgers, bins, powers, SELECTIONS, args.output_dir, ROOT, args.input, provenance)
    draw_groups(summary, args.output_dir)
    qa = draw_qa(ledgers, summary, args.output_dir, bins)
    print(json.dumps(dict(status=summary['status'], output=str(args.output_dir), channels=len(summary['channels']),
                          matched_signal_events=qa['matched_signal_events'], signal_events=qa['signal_events'])))


if __name__ == '__main__':
    main()
