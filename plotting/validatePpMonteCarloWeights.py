#!/usr/bin/env python3
"""Check ppRef MC vertex reweighting against the supplied processed pp data.

Uses the same one-dimensional, unit-area procedure as the OO validation.
This derives and checks a candidate on the same input, not independent closure.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

from validateOOMonteCarloWeights import (
    ROOT, REPOSITORY, fit_ratio, metrics, normalized, np, plt,
    read_card, read_histogram,
)
from numpy.polynomial import Polynomial


def source_parameters(source):
    block = source.read_text().split("else if(fIsPp2024Data)", 1)[1]
    match = re.search(r"fVzWeightFunction->SetParameters\(([^;]+)\);", block)
    if not match:
        raise ValueError("Cannot locate ppRef vertex parameters")
    return [float(value.strip()) for value in match.group(1).split(",")]


def run(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = REPOSITORY / "src/EECAnalyzer.cxx"
    parameters = source_parameters(source)
    inputs, histograms, weighted = {}, {}, {}
    for label, path in [("data", args.data), ("mc", args.mc)]:
        root_file = ROOT.TFile.Open(str(path), "READ")
        if not root_file or root_file.IsZombie():
            raise ValueError(f"Cannot open {path}")
        inputs[label] = {"path": str(path.resolve()),
                         "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "card": read_card(root_file)}
        histograms[label] = read_histogram(root_file, "vertexZ")
        weighted[label] = read_histogram(root_file, "vertexZweighted")
        root_file.Close()

    # Data and MC use different trigger selections; MC has no trigger requirement.
    if inputs["data"]["card"].get("ZVertexCut") != inputs["mc"]["card"].get("ZVertexCut"):
        raise ValueError("Incompatible data/MC card setting: ZVertexCut")
    cut = inputs["mc"]["card"].get("ZVertexCut", [])
    if cut != [15.0]:
        raise ValueError(f"Expected |vz| < 15 cm inputs, got ZVertexCut={cut}")
    mh, dh = histograms["mc"], weighted["data"]
    for histogram in [dh, weighted["mc"], histograms["data"]]:
        if not np.array_equal(mh["edges"], histogram["edges"]):
            raise ValueError("Incompatible vertex histogram binning")
    lower, upper = -15.0, 15.0
    centers = (mh["edges"][:-1] + mh["edges"][1:]) / 2
    select = (centers >= lower) & (centers < upper)
    x = centers[select]
    m, d = mh["counts"][select], dh["counts"][select]
    if np.any(m <= 0) or np.any(d <= 0):
        raise ValueError("Empty vertex bins; do not silently discard unsupported regions")
    mn, dn = normalized(m), normalized(d)
    ratio = dn / mn
    ratio_error = ratio * np.sqrt((dh["errors"][select] / d)**2 + (mh["errors"][select] / m)**2)
    current = Polynomial(parameters)
    refit, fit_report = fit_ratio(x, ratio, ratio_error, lower, upper)
    current_weight, candidate_weight = current(x), refit(x)
    current_prediction = normalized(m * current_weight)
    candidate_prediction = normalized(m * candidate_weight)
    supplied_prediction = normalized(weighted["mc"]["counts"][select])
    dense = current(np.linspace(lower, upper, 10001))
    comparison = {
        "range": [lower, upper], "mc_sum_in_range": float(m.sum()), "data_sum_in_range": float(d.sum()),
        "mc_entries_all_bins": mh["entries"], "data_entries_all_bins": dh["entries"],
        "raw": metrics(mn, dn), "as_supplied_weighted": metrics(supplied_prediction, dn),
        "current_function_prediction": metrics(current_prediction, dn),
        "candidate_refit_prediction": metrics(candidate_prediction, dn),
        "current_minimum_weight_on_dense_grid": float(dense.min()),
        "current_maximum_weight_on_dense_grid": float(dense.max()),
        "candidate_refit": fit_report,
        "max_supplied_raw_weighted_bin_difference": {
            label: float(np.max(np.abs(weighted[label]["counts"] - histograms[label]["counts"])))
            for label in ["data", "mc"]},
    }
    report = {
        "inputs": inputs, "analyzer_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "current_parameters_constant_first": parameters,
        "method": "Unit-area data/MC vertex comparison within |vz| < 15 cm, evaluated at 0.5 cm bin centers; error-weighted pol6 fit to normalized data/MC, using the same bins to derive and check the candidate.",
        "limitations": [
            "The candidate uses the same sample for derivation and checking; this is not independent closure.",
            "The MC input is GenGen while the data input is RecoReco; reconstruction-mode dependence requires a corresponding reconstructed MC check.",
            "Generator-event/cross-section normalization and event-level combinations of weights cannot be validated from processed one-dimensional histograms.",
            "Bin-center application approximates the unknown distribution within each 0.5 cm bin.",
            "Fit chi2 uses diagonal ratio errors; normalization correlations and fitted-weight uncertainties are not included.",
            "Both supplied weighted vertex histograms equal their raw counterparts, so this compares offline predictions rather than a weighted analyzer rerun.",
        ],
        "distributions": {"vertexZ": comparison},
    }
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                    gridspec_kw={"height_ratios": [2, 1]})
    top.errorbar(x, dn, yerr=dh["errors"][select] / d.sum(), fmt="k.", markersize=3, label="pp2024 data")
    for values, label, color in [(mn, "Supplied MC (unit weights)", "#999999"),
                                 (current_prediction, "Existing pol6, marginal check", "#c44e52"),
                                 (candidate_prediction, "Refit, same marginal input", "#0072b2")]:
        top.plot(x, values, color=color, label=label)
        bottom.plot(x, values / dn, color=color)
    top.set_ylabel("Fraction per 0.5 cm bin")
    top.legend(fontsize=9)
    top.set_title("ppRef vertex-weight check: one-dimensional, same-sample reweighting")
    bottom.axhline(1, color="black", linewidth=0.8)
    bottom.set_ylabel("MC / data")
    bottom.set_xlabel(r"$v_z$ [cm]")
    bottom.grid(alpha=0.25)
    fig.tight_layout()
    for extension in ["pdf", "png"]:
        fig.savefig(args.output_dir / f"ppEventWeightValidation.{extension}", dpi=180)
    plt.close(fig)
    (args.output_dir / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    rows = [{"bin_center": float(x[i]), "data_fraction": float(dn[i]), "raw_mc_fraction": float(mn[i]),
             "target_data_over_mc": float(ratio[i]), "target_ratio_diagonal_error": float(ratio_error[i]),
             "current_weight": float(current_weight[i]), "current_mc_over_data": float(current_prediction[i] / dn[i]),
             "candidate_weight": float(candidate_weight[i]), "candidate_mc_over_data": float(candidate_prediction[i] / dn[i])}
            for i in range(len(x))]
    with (args.output_dir / "binComparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(comparison, indent=2))
    print(f"Saved validation to {args.output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mc", type=Path, default=REPOSITORY / "merged_eecppMC2024_GenGen_processed.root")
    parser.add_argument("--data", type=Path, default=REPOSITORY / "merged_eecpp2024_run3_v7_processed.root")
    parser.add_argument("--output-dir", type=Path, default=REPOSITORY / "validation/ppEventWeights")
    run(parser.parse_args())
