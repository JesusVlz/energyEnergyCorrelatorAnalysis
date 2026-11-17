#!/usr/bin/env python3
"""Plot the applied OO veto support, validating a merged snapshot against its payload.

Run with figures/OO_ppRef_background_20261108/run_in_cmssw.sh python3
plotting/plotOOJetVetoMap.py. Nonzero payload bins veto, regardless of the
amplitude that hadd assigns to a repeated configured-map snapshot.
"""
import argparse
from array import array
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def edges(axis):
    return [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)]


def values(hist):
    return [hist.GetBinContent(x, y)
            for x in range(hist.GetNbinsX() + 2)
            for y in range(hist.GetNbinsY() + 2)]


def map_digest(hist):
    encoded = json.dumps({'x_edges': edges(hist.GetXaxis()),
                          'y_edges': edges(hist.GetYaxis()),
                          'values_x_major_including_flow': values(hist)},
                         sort_keys=True, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(encoded.encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--provenance', type=Path,
                        default=REPO / 'figures/OO_MC_firstReview/jetQA/style_provenance.json')
    parser.add_argument('--output-prefix', type=Path,
                        default=REPO / 'figures/OO_MC_firstReview/jetQA/OO_MC_jet_veto_map')
    args = parser.parse_args()
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaperSize(20, 20)
    provenance = json.loads(args.provenance.read_text())
    source_path = Path(provenance['inputs']['integrated']['path'])
    source_hash = sha256(source_path)
    if source_hash != provenance['inputs']['integrated']['sha256']:
        raise RuntimeError('JetQA input differs from its recorded provenance')
    source_file = ROOT.TFile.Open(str(source_path), 'READ')
    if not source_file or source_file.IsZombie():
        raise RuntimeError('Cannot open jetQA source')
    source_key = 'run3JetQA/configuredJetVetoMap'
    source = source_file.Get(source_key)
    if not source or not source.InheritsFrom('TH2'):
        raise RuntimeError('Configured two-dimensional veto-map snapshot is absent')
    payload_path = Path(source_file.Get('run3JetQA/VetoMapFile').GetTitle())
    if not payload_path.is_absolute():
        payload_path = REPO / payload_path
    payload_key = source_file.Get('run3JetQA/VetoMapKey').GetTitle()
    payload_file = ROOT.TFile.Open(str(payload_path), 'READ')
    if not payload_file or payload_file.IsZombie():
        raise RuntimeError('Cannot open the pinned veto payload')
    payload = payload_file.Get(payload_key)
    if not payload or not payload.InheritsFrom('TH2'):
        raise RuntimeError('Pinned two-dimensional veto payload is absent')
    if edges(source.GetXaxis()) != edges(payload.GetXaxis()) or edges(source.GetYaxis()) != edges(payload.GetYaxis()):
        raise RuntimeError('Configured and pinned map axes differ')
    source_values, payload_values = values(source), values(payload)
    if any(not math.isfinite(v) for v in source_values + payload_values):
        raise RuntimeError('Nonfinite veto-map value')
    if any((a != 0) != (b != 0) for a, b in zip(source_values, payload_values)):
        raise RuntimeError('Configured and pinned veto supports differ')
    scale_factors = sorted({a / b for a, b in zip(source_values, payload_values) if b != 0})

    mask = source.Clone('veto_mask')
    mask.SetDirectory(0)
    mask.Reset('ICES')
    vetoed = 0
    for x in range(source.GetNbinsX() + 2):
        for y in range(source.GetNbinsY() + 2):
            value = int(source.GetBinContent(x, y) != 0)
            mask.SetBinContent(x, y, value)
            mask.SetBinError(x, y, 0)
            if 1 <= x <= source.GetNbinsX() and 1 <= y <= source.GetNbinsY():
                vetoed += value
    mask.SetTitle('')
    mask.SetStats(False)
    mask.SetMinimum(0)
    mask.SetMaximum(1)
    mask.SetContour(2, array('d', [0, .5]))
    red = ROOT.TColor.GetColor('#d62728')
    ROOT.gStyle.SetPalette(2, array('i', [ROOT.kWhite, red]))
    for axis in (mask.GetXaxis(), mask.GetYaxis()):
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(.048)
        axis.SetLabelSize(.037)
        axis.SetLabelOffset(.009)
        axis.CenterTitle(True)
        axis.SetNdivisions(510)
    mask.GetXaxis().SetTitle('Detector jet #eta')
    mask.GetYaxis().SetTitle('Detector jet #phi (rad)')
    mask.GetXaxis().SetTitleOffset(1.16)
    mask.GetYaxis().SetTitleOffset(1.28)
    canvas = ROOT.TCanvas('OO_jet_veto_map_canvas', 'OO jet veto map', 900, 850)
    canvas.SetCanvasSize(900, 850)
    canvas.SetLeftMargin(.14)
    canvas.SetRightMargin(.04)
    canvas.SetTopMargin(.18)
    canvas.SetBottomMargin(.13)
    canvas.SetTicks(1, 1)
    canvas.SetFrameLineWidth(2)
    mask.Draw('COL')
    labels = ROOT.TLatex()
    labels.SetNDC(True)
    labels.SetTextFont(42)
    labels.SetTextSize(.036)
    labels.DrawLatex(.14, .944, '#bf{OO 5.36 TeV}')
    labels.DrawLatex(.14, .888, 'Jet veto map')
    allowed_box = ROOT.TBox()
    allowed_box.SetFillColor(ROOT.kWhite)
    allowed_box.SetLineColor(ROOT.kBlack)
    veto_box = ROOT.TBox()
    veto_box.SetFillColor(red)
    veto_box.SetLineColor(ROOT.kBlack)
    legend = ROOT.TLegend(.65, .864, .94, .972)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(.033)
    legend.SetMargin(.17)
    legend.AddEntry(allowed_box, 'Not vetoed', 'f')
    legend.AddEntry(veto_box, 'Vetoed', 'f')
    legend.Draw()
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()

    prefix = args.output_prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    output_paths = {extension: Path(str(prefix) + '.' + extension)
                    for extension in ('png', 'pdf', 'root', 'json')}
    for extension in ('png', 'pdf'):
        canvas.SaveAs(str(output_paths[extension]))
    raw = source.Clone('configuredJetVetoMap_raw_snapshot')
    raw.SetDirectory(0)
    pinned = payload.Clone('jetvetomap_all_pinned_payload')
    pinned.SetDirectory(0)
    numerical = ROOT.TFile(str(output_paths['root']), 'RECREATE')
    raw.Write()
    pinned.Write()
    mask.Write()
    ROOT.TNamed('source_file', str(source_path)).Write()
    ROOT.TNamed('source_key', source_key).Write()
    ROOT.TNamed('source_sha256', source_hash).Write()
    ROOT.TNamed('mask_definition', '1 iff configured payload value != 0; 0 otherwise. Binary support validated against pinned original payload. This is not jet occupancy.').Write()
    canvas.Write('canvas')
    numerical.Close()

    saved = ROOT.TFile.Open(str(output_paths['root']), 'READ')
    saved_mask = saved.Get('veto_mask')
    saved_raw = saved.Get('configuredJetVetoMap_raw_snapshot')
    if values(saved_mask) != [float(v != 0) for v in source_values]:
        raise RuntimeError('Saved mask differs from configured veto support')
    if values(saved_raw) != source_values:
        raise RuntimeError('Saved raw snapshot differs from input')
    if map_digest(saved_mask) != map_digest(mask):
        raise RuntimeError('Saved mask axes or values changed')
    saved.Close()
    summary = {
        'status': 'COMPLETE',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'script': str(Path(__file__).resolve()),
        'script_sha256': sha256(__file__),
        'source': {'path': str(source_path), 'sha256': source_hash, 'key': source_key,
                   'histogram_sha256': map_digest(source)},
        'pinned_payload': {'path': str(payload_path), 'sha256': sha256(payload_path),
                           'key': payload_key, 'histogram_sha256': map_digest(payload)},
        'mask_definition': 'Vetoed iff source bin content != 0, matching Run3JetSelection.cxx IsVetoed; map is independent of event centrality.',
        'merged_snapshot_note': 'hadd sums repeated snapshot amplitudes. Only nonzero support is plotted, after exact support validation against the original pinned payload.',
        'snapshot_over_payload_nonzero_scale_factors': scale_factors,
        'bins': {'eta': source.GetNbinsX(), 'phi': source.GetNbinsY(),
                 'regular': source.GetNbinsX() * source.GetNbinsY(), 'vetoed': vetoed,
                 'checked_including_flow': len(source_values)},
        'eta_range': [source.GetXaxis().GetXmin(), source.GetXaxis().GetXmax()],
        'phi_range_rad': [source.GetYaxis().GetXmin(), source.GetYaxis().GetXmax()],
        'source_regular_values': dict(Counter(source.GetBinContent(x, y)
                                            for x in range(1, source.GetNbinsX() + 1)
                                            for y in range(1, source.GetNbinsY() + 1))),
        'checks': {'source_file_hash_matches_jetqa_provenance': True,
                   'snapshot_support_matches_pinned_payload': True,
                   'saved_mask_values_and_axes_exact': True,
                   'saved_raw_snapshot_values_exact': True},
        'outputs': {ext: {'path': str(path), 'sha256': sha256(path)}
                    for ext, path in output_paths.items() if ext != 'json'},
    }
    output_paths['json'].write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'status': 'COMPLETE', 'output_prefix': str(prefix),
                      'checked_bins_including_flow': len(source_values),
                      'vetoed_bins': vetoed, 'snapshot_scale_factors': scale_factors}))


if __name__ == '__main__':
    main()
