#!/usr/bin/env python3
"""Check OO event-weight marginals using processed MC and data ROOT files.

Requires PyROOT, numpy and matplotlib. This is a same-sample, one-dimensional
reweighting check, not event-level closure of w(vz)*w(centrality)*w(generator).
The input ROOT files and the analyzer coefficients are never modified.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial
import ROOT


ROOT.gROOT.SetBatch(True)
REPOSITORY = Path(__file__).resolve().parents[1]


def read_histogram(root_file, name):
    histogram = root_file.Get(name)
    if not histogram or not histogram.InheritsFrom("TH1") or histogram.GetDimension() != 1:
        raise ValueError(f"Missing one-dimensional histogram: {name}")
    n = histogram.GetNbinsX()
    return {
        "edges": np.array([histogram.GetBinLowEdge(i) for i in range(1, n + 2)]),
        "counts": np.array([histogram.GetBinContent(i) for i in range(1, n + 1)]),
        "errors": np.array([histogram.GetBinError(i) for i in range(1, n + 1)]),
        "entries": float(histogram.GetEntries()),
        "underflow": float(histogram.GetBinContent(0)),
        "overflow": float(histogram.GetBinContent(n + 1)),
    }


def read_card(root_file):
    directory = root_file.Get("JCard")
    if not directory:
        return {}
    card = {}
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TObjString"):
            card[key.GetName()] = str(obj.GetString())
        elif obj.ClassName().startswith("TVectorT<"):
            # JCard vectors have a lower index of 1; Python [] adds a separate
            # zero-based bounds check. Call ROOT's operator() explicitly.
            card[key.GetName()] = [float(obj(i)) for i in range(obj.GetLwb(), obj.GetUpb() + 1)]
    return card


def coefficients(source):
    text = source.read_text()
    oo_block = text.split("} else if(fIsOOData){", 1)[1].split("else if(fIsPp2024Data)", 1)[0]
    result = {}
    for label, variable in [("vertexZ", "fVzWeightFunction"),
                            ("central", "fCentralityWeightFunctionCentral"),
                            ("peripheral", "fCentralityWeightFunctionPeripheral")]:
        match = re.search(re.escape(variable) + r"->SetParameters\(([^;]+)\);", oo_block)
        if not match:
            raise ValueError(f"Cannot locate OO parameters for {variable}")
        result[label] = [float(x.strip()) for x in match.group(1).split(",")]
    return result


def metrics(candidate, target):
    residual = candidate / target - 1
    return {
        "data_weighted_relative_rms_percent": float(100 * np.sqrt(np.sum(target * residual**2))),
        "maximum_absolute_relative_residual_percent": float(100 * np.max(np.abs(residual))),
        "total_variation_percent": float(50 * np.sum(np.abs(candidate - target))),
    }


def normalized(values):
    if not np.all(np.isfinite(values)) or np.any(values < 0) or values.sum() <= 0:
        raise ValueError("Nonfinite, negative, or empty histogram/weighted prediction")
    return values / values.sum()


def fit_ratio(x, ratio, error, lower, upper):
    select = (x >= lower) & (x < upper)
    if np.count_nonzero(select) <= 7 or np.any(error[select] <= 0):
        raise ValueError("Insufficient populated bins or invalid ratio errors for pol6 fit")
    # Fit in scaled coordinates to avoid the poorly conditioned raw c**6 basis.
    fit = Polynomial.fit(x[select], ratio[select], 6, w=1 / error[select], domain=[lower, upper])
    parameters = fit.convert().coef.tolist()
    dense = fit(np.linspace(lower, upper, 10001))
    residual = (fit(x[select]) - ratio[select]) / error[select]
    return fit, {
        "parameters_constant_first": parameters,
        "range": [lower, upper],
        "minimum_weight_on_dense_grid": float(dense.min()),
        "maximum_weight_on_dense_grid": float(dense.max()),
        "diagonal_error_chi2": float(np.sum(residual**2)),
        "ndof": int(np.count_nonzero(select) - 7),
    }


def fit_centrality(x, ratio, error):
    """Simultaneously fit two pol6 functions with a common value at 30%."""
    matrix = np.zeros((len(x), 13))
    matrix[:, 0] = 1
    for select, lower, upper, endpoint, offset in [
        (x < 30, 0, 30, 1, 1), (x >= 30, 30, 90, -1, 7),
    ]:
        scaled = 2 * (x[select] - lower) / (upper - lower) - 1
        for power in range(1, 7):
            matrix[select, offset + power - 1] = scaled**power - endpoint**power
    parameters = np.linalg.lstsq(matrix / error[:, None], ratio / error, rcond=None)[0]
    functions, details = [], {}
    for label, lower, upper, endpoint, offset in [
        ("central", 0, 30, 1, 1), ("peripheral", 30, 90, -1, 7),
    ]:
        scaled_x = Polynomial([-1 - 2 * lower / (upper - lower), 2 / (upper - lower)])
        function = Polynomial([parameters[0]])
        for power in range(1, 7):
            function += parameters[offset + power - 1] * (scaled_x**power - endpoint**power)
        dense = function(np.linspace(lower, upper, 10001))
        functions.append(function)
        details[label] = {"parameters_constant_first": function.coef.tolist(), "range": [lower, upper],
                          "minimum_weight_on_dense_grid": float(dense.min()),
                          "maximum_weight_on_dense_grid": float(dense.max())}
    details["constraint"] = "Both pol6 functions share the same value at centrality 30%; slopes are unconstrained."
    details["diagonal_error_chi2"] = float(np.sum(((matrix @ parameters - ratio) / error)**2))
    details["ndof"] = len(x) - 13
    return functions[0], functions[1], details


def run(args):
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    source = REPOSITORY / "src/EECAnalyzer.cxx"
    current = coefficients(source)
    files = {}
    inputs = {}
    for label, path in [("mc", args.mc), ("data", args.data)]:
        f = ROOT.TFile.Open(str(path), "READ")
        if not f or f.IsZombie():
            raise ValueError(f"Cannot open {path}")
        files[label] = f
        inputs[label] = {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "card": read_card(f)}

    report = {
        "inputs": inputs,
        "analyzer_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "current_parameters_constant_first": current,
        "method": "Unit-area data/MC comparison in each stated range; fixed binwise application of the current functions; pol6 refits to the same marginal inputs, with continuity enforced between the centrality pieces.",
        "limitations": [
            "No event-level vz-by-centrality distribution is available: the product of both corrections cannot be validated from these marginals.",
            "Generator-event/cross-section weights cannot be recovered or validated from the processed one-dimensional histograms.",
            "The candidate refit uses the same sample for derivation and checking; it is not independent closure.",
            "Vertex weights are evaluated at bin centers, approximating the unknown distribution within each 0.5 cm bin.",
            "Fit chi2 uses diagonal ratio errors; normalization correlations and fitted-weight uncertainties are not included.",
        ],
        "distributions": {},
    }
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex="col", gridspec_kw={"height_ratios": [2, 1]})
    rows = []
    for column, (name, weighted_name, lower, upper, xlabel) in enumerate([
        ("vertexZ", "vertexZweighted", -15, 15, r"$v_z$ [cm]"),
        ("centrality", "centralityWeighted", 0, 90, "Centrality [%]"),
    ]):
        mh, dh = read_histogram(files["mc"], name), read_histogram(files["data"], weighted_name)
        supplied = read_histogram(files["mc"], weighted_name)
        if not np.array_equal(mh["edges"], dh["edges"]) or not np.array_equal(mh["edges"], supplied["edges"]):
            raise ValueError(f"Incompatible {name} binning")
        centers = (mh["edges"][:-1] + mh["edges"][1:]) / 2
        select = (centers >= lower) & (centers < upper)
        x_plot = centers[select]
        # ForestReader::GetCentrality() = hiBin/2.0; for a 0.5%-wide
        # histogram each discrete event coordinate is the bin's low edge.
        x = mh["edges"][:-1][select] if name == "centrality" else x_plot
        if name == "centrality" and not np.allclose(np.diff(mh["edges"]), 0.5):
            raise ValueError("Centrality reconstruction requires the framework's 0.5%-wide bins")
        m, d = mh["counts"][select], dh["counts"][select]
        if np.any(m <= 0) or np.any(d <= 0):
            raise ValueError(f"Empty {name} bins: do not silently drop unsupported regions")
        mn, dn = normalized(m), normalized(d)
        ratio = dn / mn
        ratio_error = ratio * np.sqrt((dh["errors"][select] / d)**2 + (mh["errors"][select] / m)**2)
        if name == "vertexZ":
            weight = Polynomial(current[name])(x)
            refit, fit_report = fit_ratio(x, ratio, ratio_error, lower, upper)
            fitted_weight = refit(x)
            fit_report = {"vertex": fit_report}
        else:
            weight = np.where(x < 30, Polynomial(current["central"])(x), Polynomial(current["peripheral"])(x))
            central, peripheral, fit_report = fit_centrality(x, ratio, ratio_error)
            fitted_weight = np.where(x < 30, central(x), peripheral(x))
            fit_report["current_junction_30_percent"] = [float(Polynomial(current[key])(30)) for key in ["central", "peripheral"]]
            fit_report["candidate_junction_30_percent"] = [float(central(30)), float(peripheral(30))]
        predicted = normalized(m * weight)
        fitted = normalized(m * fitted_weight)
        as_supplied = normalized(supplied["counts"][select])
        detail = {
            "range": [lower, upper], "mc_sum_in_range": float(m.sum()), "data_sum_in_range": float(d.sum()),
            "mc_entries_all_bins": mh["entries"], "data_entries_all_bins": dh["entries"],
            "max_supplied_raw_weighted_bin_difference": float(np.max(np.abs(supplied["counts"] - mh["counts"]))),
            "raw": metrics(mn, dn), "as_supplied_weighted": metrics(as_supplied, dn),
            "current_function_prediction": metrics(predicted, dn), "candidate_refit_prediction": metrics(fitted, dn),
            "current_weight_range_at_event_coordinates": [float(weight.min()), float(weight.max())],
            "candidate_refit": fit_report,
        }
        if name == "centrality":
            detail["region_fractions_percent"] = []
            for a, b in [(0, 30), (30, 50), (50, 90)]:
                region = (x >= a) & (x < b)
                detail["region_fractions_percent"].append({"range": [a, b], **{
                    label: float(100 * values[region].sum()) for label, values in
                    [("data", dn), ("raw_mc", mn), ("current_prediction", predicted), ("candidate_prediction", fitted)]}})
        report["distributions"][name] = detail
        top, bottom = axes[0, column], axes[1, column]
        top.plot(x_plot, dn, "k.", label="OO data")
        for values, label, color in [(mn, "Supplied MC (unit weights)", "#999999"),
                                      (predicted, "Current function, marginal check", "#c44e52"),
                                      (fitted, "Refit, same marginal input", "#0072b2")]:
            top.plot(x_plot, values, color=color, label=label)
            bottom.plot(x_plot, values / dn, color=color)
        top.set_ylabel("Fraction per bin")
        top.legend(fontsize=8)
        top.set_title("Vertex weight" if name == "vertexZ" else "Centrality weight (event coordinates = bin low edges)")
        bottom.axhline(1, color="black", linewidth=0.8)
        bottom.set_ylabel("MC / data")
        bottom.set_xlabel(xlabel)
        bottom.grid(alpha=0.25)
        for i in range(len(x)):
            rows.append({"distribution": name, "event_coordinate": float(x[i]), "bin_center": float(x_plot[i]),
                         "data_fraction": float(dn[i]), "raw_mc_fraction": float(mn[i]),
                         "target_data_over_mc": float(ratio[i]), "target_ratio_diagonal_error": float(ratio_error[i]),
                         "current_weight": float(weight[i]), "current_mc_over_data": float(predicted[i] / dn[i]),
                         "candidate_weight": float(fitted_weight[i]), "candidate_mc_over_data": float(fitted[i] / dn[i])})
        print(name, json.dumps({key: detail[key] for key in ["raw", "current_function_prediction", "candidate_refit_prediction"]}))
    fig.suptitle("OO event-weight check: one-dimensional, same-sample reweighting")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for extension in ["pdf", "png"]:
        fig.savefig(output / f"ooEventWeightValidation.{extension}", dpi=180)
    plt.close(fig)
    (output / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (output / "binComparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for root_file in files.values():
        root_file.Close()
    print(f"Saved validation to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mc", type=Path, default=REPOSITORY / "eecAnalysis_OOMC_eventWeights_merged_processed.root")
    parser.add_argument("--data", type=Path, default=REPOSITORY / "merged_eecOO2025_processed_121026.root")
    parser.add_argument("--output-dir", type=Path, default=REPOSITORY / "validation/ooEventWeights")
    run(parser.parse_args())
