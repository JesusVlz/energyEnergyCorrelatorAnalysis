#!/usr/bin/env python3
"""Inspect local ROOT candidates without installing or selecting any payload.

Requires PyROOT. Optional plots compare histogram contents, not physics closure.
"""
import argparse
import collections
import hashlib
import json
import math
from pathlib import Path


def inspect(path, ROOT, plots_dir=None):
    path = Path(path)
    result = dict(filename=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  size_bytes=path.stat().st_size, root_version=ROOT.gROOT.GetVersion(),
                  keys=[], axis_definitions={})
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie():
        raise ValueError(f'Unreadable ROOT file: {path}')
    try:
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            record = dict(name=key.GetName(), cycle=key.GetCycle(),
                          **{'class': key.GetClassName()}, title=obj.GetTitle())
            if obj.InheritsFrom('TH1'):
                record.update(dimension=obj.GetDimension(), entries=obj.GetEntries(), axes={})
                axes = [('x', obj.GetXaxis()), ('y', obj.GetYaxis()), ('z', obj.GetZaxis())]
                for name, axis in axes[:obj.GetDimension()]:
                    definition = dict(name=axis.GetName(), title=axis.GetTitle(),
                                      bins=axis.GetNbins(), min=axis.GetXmin(), max=axis.GetXmax(),
                                      edges=[axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2)])
                    label = next((k for k, v in result['axis_definitions'].items() if v == definition), None)
                    if label is None:
                        label = f"axis_{len(result['axis_definitions'])+1}"
                        result['axis_definitions'][label] = definition
                    record['axes'][name] = label
                if obj.GetDimension() == 2:
                    regular = [obj.GetBinContent(i, j) for i in range(1, obj.GetNbinsX()+1)
                               for j in range(1, obj.GetNbinsY()+1)]
                    flows = [obj.GetBinContent(i, j) for i in range(obj.GetNbinsX()+2)
                             for j in range(obj.GetNbinsY()+2)
                             if i in (0, obj.GetNbinsX()+1) or j in (0, obj.GetNbinsY()+1)]
                    record['regular_bins'] = dict(total=len(regular), nonzero=sum(v != 0 for v in regular),
                                                 nonfinite=sum(not math.isfinite(v) for v in regular),
                                                 min=min(regular), max=max(regular))
                    counts = collections.Counter(regular)
                    if len(counts) <= 25:
                        record['regular_bins']['value_counts'] = {str(k): v for k, v in sorted(counts.items())}
                    record['flow_bins'] = dict(total=len(flows), nonzero=sum(v != 0 for v in flows))
            result['keys'].append(record)
        default, all_regions = f.Get('jetvetomap'), f.Get('jetvetomap_all')
        if default and all_regions:
            for a, b in [(default.GetXaxis(), all_regions.GetXaxis()), (default.GetYaxis(), all_regions.GetYaxis())]:
                if a.GetNbins() != b.GetNbins() or any(a.GetBinLowEdge(i) != b.GetBinLowEdge(i)
                                                     for i in range(1, a.GetNbins()+2)):
                    raise ValueError('Cannot compare maps with different axes')
            different = [(i, j) for i in range(1, default.GetNbinsX()+1)
                         for j in range(1, default.GetNbinsY()+1)
                         if default.GetBinContent(i, j) != all_regions.GetBinContent(i, j)]
            result['key_comparison'] = dict(keys=['jetvetomap', 'jetvetomap_all'], different_regular_bins=len(different))
            if plots_dir is not None:
                plots_dir = Path(plots_dir)
                plots_dir.mkdir(parents=True, exist_ok=True)
                canvas = ROOT.TCanvas('audit_'+path.stem, 'Candidate content audit', 1500, 500)
                canvas.Divide(3, 1)
                difference = all_regions.Clone('key_difference')
                difference.SetDirectory(0)
                difference.Reset()
                for i, j in different:
                    difference.SetBinContent(i, j, 1)
                for panel, hist, title in [
                        (1, all_regions, 'jetvetomap_all (user/HIN helper key)'),
                        (2, default, 'jetvetomap (pp correctionlib key)'),
                        (3, difference, 'Keys disagree (1 = different bin)')]:
                    canvas.cd(panel)
                    ROOT.gPad.SetRightMargin(0.16)
                    ROOT.gPad.SetLeftMargin(0.13)
                    ROOT.gPad.SetBottomMargin(0.14)
                    hist.SetStats(False)
                    hist.SetMinimum(0)
                    hist.SetMaximum(1 if panel == 3 else 100)
                    for axis in (hist.GetXaxis(), hist.GetYaxis()):
                        axis.SetTitleFont(42)
                        axis.SetLabelFont(42)
                        axis.SetTitleSize(0.045)
                        axis.SetTitleOffset(1.2)
                    hist.SetTitle(title+';jet #eta;jet #phi [rad]')
                    hist.Draw('COLZ')
                plot = plots_dir/(path.stem+'_key_comparison.png')
                canvas.SaveAs(str(plot))
                result['validation_plots'] = [dict(path=str(plot), sha256=hashlib.sha256(plot.read_bytes()).hexdigest(),
                                                  scope='Histogram contents only; no event/physics validation')]
    finally:
        f.Close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('maps', nargs='+', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--plots-dir', type=Path)
    args = parser.parse_args()
    if args.output and args.output.resolve() in [p.resolve() for p in args.maps]:
        parser.error('Metadata output cannot overwrite an input map')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    data = json.dumps([inspect(p, ROOT, args.plots_dir) for p in args.maps], indent=2, allow_nan=False)+'\n'
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end='')


if __name__ == '__main__':
    main()
