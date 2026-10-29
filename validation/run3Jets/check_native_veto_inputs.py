#!/usr/bin/env python3
"""Replay the C++ HIN Jet ID and veto on retained, bounded native jet ledgers.

This is unweighted jet-kinematics validation. It does not run the analyzer,
apply its event/trigger selection, establish complete >15 GeV jet retention,
validate a correction chain, or establish full-campaign map applicability.
Both rawpt and stored jtpt are evaluated; neither is chosen as nominal.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import shlex
import subprocess
import tempfile
from datetime import datetime, timezone


SAMPLES = (
    "OO_data_production", "OO_MC_local", "ppRef_data_production",
    "ppRef_MC_production", "ppRef_MC_local",
)
ARRAYS = (
    "rawpt", "jtpt", "jteta", "jtphi", "jtPfCHF", "jtPfNHF", "jtPfCEF",
    "jtPfNEF", "jtPfMUF", "jtPfCHM", "jtPfNHM", "jtPfCEM", "jtPfNEM", "jtPfMUM",
)
COUNTS = {"jtPfCHM", "jtPfNHM", "jtPfCEM", "jtPfNEM", "jtPfMUM"}
MAPS = {
    "OO2025": "Summer24Prompt25_RunCDEFG.root",
    "ppRef2024": "Summer24Prompt24_RunBCDEFGHI.root",
}
LIMITATIONS = [
    "Denominators are only the retained bounded JSONL entries, including entries with no retained reconstructed jets.",
    "No analyzer event cuts, trigger cuts, centrality cuts, track cuts, pThat cuts, weights, mixing, or response selection are applied.",
    "The replay visits every retained reconstructed jet in its original native order; it cannot recover jets omitted by production or prove complete pT>15 GeV retention.",
    "rawpt and stored jtpt are independent explicit sensitivity cases; neither is selected as nominal, and no new JEC is computed.",
    "Detector jteta/jtphi are used; the WTA coordinates are not substituted.",
    "OO uses pinned HINPbPb ID and ppRef uses pinned HINpp ID. Their differing cuts and HINpp five-category count definition are intentional.",
    "The selected era map is queried identically for data and MC, but these different bounded samples do not demonstrate equal event-loss fractions.",
    "Actual observed data run ranges are only ranges in the retained ledgers; MC run=1 is a simulation identifier, not data-era coverage.",
    "No claim of full-campaign run coverage, exact production applicability, or map approval follows from this bounded replay.",
    "OO_MC_local is a local retained candidate; equivalence to the current production campaign remains unconfirmed. No OO_MC_production event ledger is available here.",
]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_event(record, source):
    if record.get("original_order_preserved") is not True:
        raise ValueError(f"{source}: original native jet order is not attested")
    values = record["values"]
    nref = values["nref"]
    if type(nref) is not int or nref < 0:
        raise ValueError(f"{source}: invalid nref")
    for field in ARRAYS:
        array = values[field]
        if not isinstance(array, list) or len(array) != nref:
            raise ValueError(f"{source}: {field} length differs from nref")
        for value in array:
            if field in COUNTS:
                if type(value) is not int or value < 0:
                    raise ValueError(f"{source}: invalid {field} count")
            elif not math.isfinite(value):
                raise ValueError(f"{source}: nonfinite {field}")
    return values


def make_inputs(ROOT, values, index):
    jet = ROOT.EECJetSelection.JetInputs()
    for destination, source in (
        ("eta", "jteta"), ("phi", "jtphi"), ("rawPt", "rawpt"), ("storedPt", "jtpt"),
        ("neutralHadronFraction", "jtPfNHF"), ("neutralEmFraction", "jtPfNEF"),
        ("chargedHadronFraction", "jtPfCHF"), ("chargedEmFraction", "jtPfCEF"),
        ("muonFraction", "jtPfMUF"),
    ):
        setattr(jet, destination, values[source][index])
    jet.chargedMultiplicity = sum(values[name][index] for name in ("jtPfCHM", "jtPfCEM", "jtPfMUM"))
    jet.neutralMultiplicity = sum(values[name][index] for name in ("jtPfNHM", "jtPfNEM"))
    jet.numberOfConstituents = jet.chargedMultiplicity + jet.neutralMultiplicity
    # Mirrors the native reader. HINpp explicitly uses these category sums;
    # the generic CHS alternatives retain the separate full-HF ambiguity guard.
    jet.multiplicitiesAreLowerBounds = True
    return jet


def inspect_map(ROOT, path, key):
    source = ROOT.TFile.Open(str(path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot inspect map {path}")
    histogram = source.Get(key)
    if not histogram or not histogram.InheritsFrom("TH2"):
        raise RuntimeError(f"Cannot inspect {key} in {path}")
    result = {
        "path": str(path), "sha256": sha256(path), "histogram_key": key,
        "root_keys": [item.GetName() for item in source.GetListOfKeys()],
        "axes": {},
    }
    for name, axis in (("eta_x", histogram.GetXaxis()), ("phi_y_radians", histogram.GetYaxis())):
        result["axes"][name] = {
            "title": axis.GetTitle(), "bins": axis.GetNbins(),
            "edges": [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)],
        }
    metadata_path = path.with_suffix(".metadata.json")
    if metadata_path.is_file():
        result["metadata_path"] = str(metadata_path)
        result["metadata_sha256"] = sha256(metadata_path)
    source.Close()
    return result


def new_counts():
    return dict.fromkeys((
        "events", "events_with_native_jets", "native_jets", "jet_id_pass_all_pt",
        "jet_id_fail_all_pt", "jets_pt_gt15", "jets_pt_gt15_and_id", "eligible_jets",
        "eligible_jets_in_veto_region", "events_with_eligible_jets", "events_rejected", "events_kept",
        "eligible_jets_after_event_veto",
    ), 0)


def make_plots(sample, cases, map_info, output, arrays):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    eta_edges = np.asarray(map_info["axes"]["eta_x"]["edges"])
    phi_edges = np.asarray(map_info["axes"]["phi_y_radians"]["edges"])
    figure, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for row, convention in enumerate(("rawpt", "jtpt")):
        points = cases[convention]
        occupancies = {}
        for category in ("before", "after", "vetoed"):
            data = np.asarray(points[category], dtype=float).reshape((-1, 2))
            occupancy, _, _ = np.histogram2d(data[:, 0], data[:, 1], bins=(eta_edges, phi_edges))
            occupancies[category] = occupancy
            arrays[f"{sample}_{convention}_{category}"] = occupancy.astype("int64")
        maximum = max(1, occupancies["before"].max())
        for column, (category, title) in enumerate((
            ("before", "eligible jets before event veto"),
            ("after", "eligible jets in surviving events"),
            ("vetoed", "eligible jets inside veto regions"),
        )):
            axis = axes[row, column]
            plot = axis.pcolormesh(eta_edges, phi_edges, occupancies[category].T,
                                   shading="flat", vmin=0, vmax=maximum)
            axis.set(xlabel="detector jet eta", ylabel="detector jet phi [rad]",
                     title=f"{convention}: {title}")
            figure.colorbar(plot, ax=axis, label="unweighted retained jets")
    figure.suptitle(f"{sample}: bounded native replay, no analyzer event/trigger selection")
    destination = output / f"{sample}_eta_phi.png"
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--samples", nargs="+", choices=SAMPLES, default=list(SAMPLES))
    parser.add_argument("--histogram-key", default="jetvetomap_all")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError(f"Output evidence directory must be empty: {output}")

    import ROOT
    ROOT.gROOT.SetBatch(True)
    # Compile the exact current implementation, then let ROOT expose only its
    # declarations. This also works on hosts whose interpreter cannot JIT all
    # std::unique_ptr implementation symbols from the .cxx directly.
    build_directory = tempfile.TemporaryDirectory(prefix="eec-native-veto-core-")
    shared_library = Path(build_directory.name) / "libRun3JetSelectionReplay.so"
    root_flags = shlex.split(subprocess.check_output(["root-config", "--cflags", "--libs"], text=True))
    build_command = ["g++", "-shared", "-fPIC", "-O2", str(repo / "src/Run3JetSelection.cxx"),
                     *root_flags, "-o", str(shared_library)]
    subprocess.run(build_command, check=True)
    if ROOT.gSystem.Load(str(shared_library)) < 0:
        raise RuntimeError("Cannot load compiled Run3JetSelection shared library")
    ROOT.gInterpreter.AddIncludePath(str(repo / "src"))
    if not ROOT.gInterpreter.Declare('#include "Run3JetSelection.h"'):
        raise RuntimeError("Cannot load the C++ Run3JetSelection declarations")
    prescriptions = {
        "OO2025": ROOT.EECJetSelection.JetIDPrescription.HINPbPb,
        "ppRef2024": ROOT.EECJetSelection.JetIDPrescription.HINpp,
    }
    map_paths = {era: repo / "jetVetoMaps" / era / filename for era, filename in MAPS.items()}
    maps = {era: ROOT.EECJetSelection.JetVetoMap(str(path), args.histogram_key) for era, path in map_paths.items()}
    map_info = {era: inspect_map(ROOT, path, args.histogram_key) for era, path in map_paths.items()}
    code_files = [Path(__file__).resolve(), repo / "src/Run3JetSelection.cxx",
                  repo / "src/Run3JetSelection.h", repo / "src/Run3JetSelectionInputs.h"]
    report = {
        "status": "bounded_native_kinematic_replay", "created_utc": datetime.now(timezone.utc).isoformat(),
        "command_argv": sys.argv, "root_version": ROOT.gROOT.GetVersion(),
        "implementation": "Current repository C++ EECJetSelection implementation compiled into a temporary shared library; declarations exposed by ROOT.gInterpreter.Declare",
        "core_build_command": build_command, "compiled_core_sha256": sha256(shared_library),
        "code_sha256": {str(path): sha256(path) for path in code_files},
        "maps": map_info, "event_weighting": "none; every retained entry has unit weight",
        "pt_conventions": {"rawpt": "native rawpt, no new correction", "jtpt": "native stored jtpt, correction provenance not inferred"},
        "nominal_pt_convention_selected": False, "histogram_key": args.histogram_key,
        "HIN_helpers_commit": "c64293a132af388af1827044779af36e57f884d2",
        "limitations": LIMITATIONS, "samples": {},
        "unavailable_source": "OO_MC_production has no retained .events.jsonl in the supplied Goal 2 evidence",
    }
    (output / "native_replay_configuration.json").write_text(
        json.dumps({"execution_completed": False, "configuration": report}, indent=2, allow_nan=False) + "\n")
    occupancy_arrays = {}
    mask_path = output / "native_replay_event_masks.jsonl"
    with mask_path.open("w") as masks:
        for sample in args.samples:
            era = "OO2025" if sample.startswith("OO_") else "ppRef2024"
            ledger = args.evidence_dir / (sample + ".events.jsonl")
            provenance_path = args.evidence_dir / (sample + ".json")
            provenance = json.loads(provenance_path.read_text())
            records = [json.loads(line) for line in ledger.read_text().splitlines() if line.strip()]
            expected = provenance["events_inspected"]
            if len(records) != expected:
                raise ValueError(f"{sample}: ledger count {len(records)} differs from original probe {expected}")
            native_values = [check_event(record, sample) for record in records]
            kinematic_ranges = {}
            for field in ("rawpt", "jtpt", "jteta", "jtphi"):
                flattened = [value for values in native_values for value in values[field]]
                kinematic_ranges[field] = {"count": len(flattened),
                    "minimum": min(flattened) if flattened else None,
                    "maximum": max(flattened) if flattened else None}
            fraction_endpoint_excursions = []
            for record, values in zip(records, native_values):
                for field in ("jtPfCHF", "jtPfNHF", "jtPfCEF", "jtPfNEF", "jtPfMUF"):
                    for index, value in enumerate(values[field]):
                        if value < 0 or value > 1:
                            fraction_endpoint_excursions.append({"entry": record["entry"], "jet_index": index,
                                "field": field, "exact_stored_value": value,
                                "handling": "unchanged value passed to core validation and all thresholds; no clipping"})
            counts = {name: new_counts() for name in ("rawpt", "jtpt")}
            points = {name: {category: [] for category in ("before", "after", "vetoed")} for name in counts}
            runs, lumis, entries = set(), set(), []
            differences = []
            for record in records:
                values = check_event(record, sample)
                jets = [make_inputs(ROOT, values, i) for i in range(values["nref"])]
                ids = [bool(ROOT.EECJetSelection.PassJetID(jet, prescriptions[era])) for jet in jets]
                runs.add(values["run"]); lumis.add((values["run"], values["lumi"])); entries.append(record["entry"])
                mask = {"sample": sample, "entry": record["entry"], "run": values["run"],
                        "lumi": values["lumi"], "event": values["evt"], "nref": len(jets), "cases": {}}
                for convention, counter in counts.items():
                    pts = values[convention]
                    eligible = [bool(ROOT.EECJetSelection.IsEligibleVetoJet(jet, pts[i], prescriptions[era])) for i, jet in enumerate(jets)]
                    vetoed = [bool(ROOT.EECJetSelection.RejectEventForJet(jet, pts[i], maps[era], prescriptions[era])) for i, jet in enumerate(jets)]
                    reject = any(vetoed)
                    counter["events"] += 1
                    counter["events_with_native_jets"] += bool(jets)
                    counter["native_jets"] += len(jets)
                    counter["jet_id_pass_all_pt"] += sum(ids)
                    counter["jet_id_fail_all_pt"] += len(jets) - sum(ids)
                    counter["jets_pt_gt15"] += sum(pt > 15 for pt in pts)
                    counter["jets_pt_gt15_and_id"] += sum(pt > 15 and passed for pt, passed in zip(pts, ids))
                    counter["eligible_jets"] += sum(eligible)
                    counter["eligible_jets_in_veto_region"] += sum(vetoed)
                    counter["events_with_eligible_jets"] += any(eligible)
                    counter["events_rejected"] += reject
                    counter["events_kept"] += not reject
                    counter["eligible_jets_after_event_veto"] += 0 if reject else sum(eligible)
                    for i, jet in enumerate(jets):
                        coordinate = (float(jet.eta), float(ROOT.EECJetSelection.DetectorPhi(jet.phi)))
                        if eligible[i]:
                            points[convention]["before"].append(coordinate)
                            if not reject: points[convention]["after"].append(coordinate)
                        if vetoed[i]: points[convention]["vetoed"].append(coordinate)
                    mask["cases"][convention] = {"reject_event": reject,
                        "jet_id_pass_indices": [i for i, flag in enumerate(ids) if flag],
                        "eligible_indices": [i for i, flag in enumerate(eligible) if flag],
                        "vetoed_indices": [i for i, flag in enumerate(vetoed) if flag]}
                if mask["cases"]["rawpt"]["reject_event"] != mask["cases"]["jtpt"]["reject_event"]:
                    differences.append(record["entry"])
                masks.write(json.dumps(mask, allow_nan=False) + "\n")
            for counter in counts.values():
                counter["bounded_retained_entry_loss_fraction"] = counter["events_rejected"] / counter["events"] if counter["events"] else None
            sample_report = {
                "era_map": era, "jet_id": "HINPbPb" if era == "OO2025" else "HINpp",
                "jet_tree": provenance["expected_tree"], "radius": 0.4,
                "input_root_file": provenance["input"], "origin": provenance["origin"],
                "production_configuration_match": provenance["production_configuration_match"],
                "input_ledger": str(ledger.resolve()), "input_ledger_sha256": sha256(ledger),
                "probe_summary": str(provenance_path.resolve()), "probe_summary_sha256": sha256(provenance_path),
                "observed_runs": sorted(runs), "observed_run_min": min(runs) if runs else None,
                "observed_run_max": max(runs) if runs else None, "observed_run_lumi_pairs": len(lumis),
                "run_identifier_is_data_coverage": "_MC_" not in sample,
                "retained_entry_min": min(entries) if entries else None, "retained_entry_max": max(entries) if entries else None,
                "retained_native_kinematic_ranges": kinematic_ranges,
                "fraction_endpoint_excursions": fraction_endpoint_excursions,
                "counts": counts, "raw_vs_stored_event_mask_difference_entries": differences,
            }
            if not args.no_plots:
                sample_report["eta_phi_plot"] = make_plots(sample, points, map_info[era], output, occupancy_arrays)
            report["samples"][sample] = sample_report
    if occupancy_arrays:
        import numpy as np
        np.savez_compressed(output / "native_replay_occupancy.npz", **occupancy_arrays)
    report["artifacts"] = {path.name: {"sha256": sha256(path), "size_bytes": path.stat().st_size}
                           for path in sorted(output.iterdir()) if path.is_file()}
    report["execution_completed"] = True
    destination = output / "native_replay_summary.json"
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"summary": str(destination), "samples": {
        sample: {convention: {key: counts[key] for key in ("events", "native_jets", "jet_id_pass_all_pt", "eligible_jets", "events_rejected")}
                 for convention, counts in info["counts"].items()} for sample, info in report["samples"].items()}}, indent=2))


if __name__ == "__main__":
    main()
