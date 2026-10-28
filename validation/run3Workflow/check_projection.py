#!/usr/bin/env python3
"""Synthetic raw inclusive projection test; no real detector inputs."""
import argparse
from array import array
from pathlib import Path
import ROOT


def vector(name, values):
    v = ROOT.TVectorT("float")(len(values))
    for i, value in enumerate(values):
        v[i] = value
    v.Write(name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scratch", type=Path)
    args = parser.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_name("projectInclusiveRaw.C").resolve()
    ROOT.gROOT.SetBatch(True)
    if ROOT.gROOT.LoadMacro(str(source)) != 0:
        raise RuntimeError("Cannot compile/load inclusive projection macro")
    for radius in (.4, .8):
        code = round(radius*10)
        first = args.scratch / f"projection_fixture_R{code:02d}.root"
        final = args.scratch / f"projection_result_R{code:02d}.root"
        if first.exists() or final.exists():
            raise RuntimeError("Use a fresh scratch directory")
        f = ROOT.TFile(str(first), "RECREATE")
        f.mkdir("JCard").cd()
        vector("InclusiveMode", [1])
        vector("WeightExponent", [1, 2])
        vector("JetRadius", [radius])
        f.mkdir("inclusive").cd()
        raw = ROOT.THnSparseD("inclusiveRawEEC", "synthetic", 6,
                              array("i", [8, 3, 2, 13, 3, 2]),
                              array("d", [0, 40, .5, -.5, -.5, -.5]),
                              array("d", [2*radius, 160, 500, 12.5, 2.5, 1.5]))
        raw.Sumw2()
        raw.SetBinEdges(1, array("d", [40, 80, 120, 160]))
        raw.SetBinEdges(2, array("d", [.5, 2, 500]))
        for j, pt in enumerate((60, 100, 140)):
            for n in (1, 2):
                # Include wide-angle pairs beyond .4. They must not be truncated.
                raw.Fill(array("d", [1.5*radius, pt, 2.5, 0, 0, n-1]), (j+1)*15**n)
        raw.Write()
        ROOT.TNamed("TrackSelection", "strict pT > 2 GeV; input threshold remains configured separately").Write()
        exposure = ROOT.TH1D("inclusiveJetExposure", "synthetic", 3, array("d", [40, 80, 120, 160]))
        exposure.Fill(60, 2)
        exposure.Write()
        f.Close()
        ROOT.projectInclusiveRaw(str(first), str(final))
        result = ROOT.TFile(str(final))
        for j, (lo, hi) in enumerate(((40, 80), (80, 120), (120, 160))):
            for n in (1, 2):
                h = result.Get(f"INC_R{code:02d}_J{lo}_{hi}_N{n}_rawSelected")
                expected = (j+1)*15**n
                if not h or abs(h.Integral()-expected) > 1e-8:
                    raise AssertionError("Projection lost raw weighted yield or nominal channel")
                if h.GetXaxis().GetXmax() < 2*radius-1e-5:
                    raise AssertionError("Projection truncated wide angular support")
        if not result.Get("correction_state") or not result.Get("JCard/JetRadius"):
            raise AssertionError("Projection lost correction/radius provenance")
        result.Close()
    print("PASS: twelve synthetic raw projections, both powers, jet bins, R04/R08 wide support and metadata; no physics closure")


if __name__ == "__main__":
    main()
