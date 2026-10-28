#!/usr/bin/env python3
"""Inspect actual capped raw outputs; never certify corrections or acceptance."""
import argparse
import json
from pathlib import Path
import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    out = {}
    for system in ("OO", "pp"):
        path = args.directory / (system+"_R04_raw.root")
        f = ROOT.TFile(str(path))
        coverage = f.Get("inclusive/inclusiveEventCoverage")
        raw = f.Get("inclusive/inclusiveRawEEC")
        exposure = f.Get("inclusive/inclusiveJetExposure")
        if not coverage or not raw or not exposure:
            raise RuntimeError("Missing real inclusive diagnostic outputs")
        stages = {coverage.GetYaxis().GetBinLabel(y):sum(coverage.GetBinContent(x,y) for x in range(1,6)) for y in range(1,5)}
        categories = {coverage.GetXaxis().GetBinLabel(x):coverage.GetBinContent(x,4) for x in range(1,6)}
        edges = f.Get("JCard/CentralityBinEdges")
        projected = ROOT.TFile(str(path.with_name(system+"_R04_projected.root")))
        spectra = {f"J{lo}_{hi}_N{n}":projected.Get(f"INC_R04_J{lo}_{hi}_N{n}_rawSelected").Integral()
                   for lo,hi in ((40,80),(80,120),(120,160)) for n in (1,2)}
        out[system] = dict(status="BOUNDED_RAW_DIAGNOSTIC_ONLY", stage_counts=stages,
                           selected_activity_categories=categories, raw_pair_entries=raw.GetEntries(),
                           jet_weighted_exposure=exposure.Integral(), projected_raw_integrals=spectra,
                           centrality_axis_edges=[edges.GetMatrixArray()[i] for i in range(edges.GetNrows())],
                           correction_state=f.Get("inclusive/CorrectionState").GetTitle(),
                           truncated_files=f.Get("inclusive/inclusiveInputStatus").GetBinContent(2))
        if stages["Read from forest"] != 100 or out[system]["truncated_files"] != 1:
            raise RuntimeError("Pilot cap/coverage incorrect")
        if raw.GetEntries() <= 0 or exposure.Integral() <= 0:
            raise RuntimeError("Empty real selected-event diagnostic")
        projected.Close()
        f.Close()
    with args.output.open("x") as stream:
        json.dump(out, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
