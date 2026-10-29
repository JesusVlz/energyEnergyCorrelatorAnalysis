#!/usr/bin/env python3
"""Check bounded MC analyzer output; this does not certify physics closure."""

import argparse
import json
import math
from pathlib import Path

from run_mc_smoke import card_values


def check_file(ROOT, args, report):
    source = ROOT.TFile.Open(str(args.root_file), "READ")
    if not source or source.IsZombie():
        raise ValueError(f"Cannot open ROOT output: {args.root_file}")
    try:
        def required(path):
            obj = source.Get(path)
            if not obj:
                raise ValueError(f"Missing ROOT object: {path}")
            return obj

        values = card_values(args.card.read_text())
        report["card_values"] = {}
        for key in ("Run3JetID", "ApplyJetVeto", "JetVetoInputMinPt", "MaxEventsPerFile", "DataType", "McCorrelationType"):
            if key not in values:
                raise ValueError(f"Missing scratch-card key: {key}")
            obj = required("JCard/" + key)
            # JCard vectors use a ROOT lower bound of 1. PyROOT's checked
            # __getitem__ treats it as a zero-based sequence and rejects [1]
            # for a single-element vector; read its first physical element.
            actual = float(obj.GetMatrixArray()[0])
            expected = float(values[key])
            if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=1e-6):
                raise ValueError(f"JCard/{key}: {actual} != card value {expected}")
            report["card_values"][key] = actual
        if any(float(values[key]) != 1 for key in ("Run3JetID", "ApplyJetVeto")):
            raise ValueError("Scratch card must enable Run3JetID and ApplyJetVeto")
        if float(values["MaxEventsPerFile"]) != args.max_events:
            raise ValueError("Scratch-card MaxEventsPerFile differs from --max-events")
        for key in ("JetVetoMapFile", "JetVetoMapHistogram", "JetVetoPtConvention"):
            obj = required("JCard/" + key + "String")
            actual = str(obj.GetString().Data())
            if actual != values.get(key):
                raise ValueError(f"Persisted {key} differs from the scratch card")
            report["card_values"][key] = actual

        # THnSparse.GetNbins is the allocated-bin count, not the dense product.
        # Scan contents and errors of stored bins, including TH1 flow bins.
        report["histograms"] = {}

        def finite_histogram(obj, path):
            entries = float(obj.GetEntries())
            if not math.isfinite(entries):
                raise ValueError(f"Nonfinite entries: {path}")
            sparse = obj.InheritsFrom("THnSparse")
            count = int(obj.GetNbins() if sparse else obj.GetNcells())
            if sparse:
                # A Python zero can also match the Int_t* coordinate overload.
                # Select linear-bin overloads explicitly to avoid null-pointer
                # queries (and recovered ROOT crashes) on allocated bin zero.
                sparse_content = obj.GetBinContent.__overload__("Long64_t, Int_t*")
                bin_content = lambda index: sparse_content(index, ROOT.nullptr)
                bin_error = obj.GetBinError.__overload__("Long64_t")
            else:
                bin_content = obj.GetBinContent
                bin_error = obj.GetBinError
            for index in range(count):
                if not math.isfinite(float(bin_content(index))) or not math.isfinite(float(bin_error(index))):
                    raise ValueError(f"Nonfinite content/error: {path}, bin {index}")
            report["histograms"][path] = {"entries": entries, "stored_bins_checked": count}

        def walk(directory, prefix=""):
            # Get(name) resolves the newest cycle; process each object once.
            for name in sorted({key.GetName() for key in directory.GetListOfKeys()}):
                obj = directory.Get(name)
                path = prefix + name
                if obj.InheritsFrom("TDirectory"):
                    walk(obj, path + "/")
                elif obj.InheritsFrom("TH1") or obj.InheritsFrom("THnSparse"):
                    finite_histogram(obj, path)

        walk(source)
        counters = required("eventQA/jetVetoEventCounts")
        if counters.GetNbinsX() != 6:
            raise ValueError("Jet-veto event counter must have six bins")
        counts = [float(counters.GetBinContent(index)) for index in range(1, 7)]
        if any(value < 0 or not value.is_integer() for value in counts):
            raise ValueError("Jet-veto event counters must be nonnegative unweighted integers")
        report["veto_counts"] = {}
        for stream, offset in (("signal", 0), ("mixing", 3)):
            checked, kept, vetoed = counts[offset:offset + 3]
            report["veto_counts"][stream] = {"checked": checked, "kept": kept, "vetoed": vetoed}
            if checked != kept + vetoed:
                raise ValueError(f"{stream}: checked != kept + vetoed")
            before = required(f"eventQA/{stream}EligibleJetEtaPhiBeforeVeto")
            after = required(f"eventQA/{stream}EligibleJetEtaPhiAfterVeto")
            if not before.InheritsFrom("TH2") or not after.InheritsFrom("TH2") or before.GetNcells() != after.GetNcells():
                raise ValueError(f"{stream}: invalid before/after eta-phi histograms")
            if any(after.GetBinContent(index) > before.GetBinContent(index) for index in range(before.GetNcells())):
                raise ValueError(f"{stream}: after-veto count exceeds before-veto count")

        events = required("nEvents")
        processed = float(events.GetBinContent(1))
        if not 0 < processed <= args.max_events:
            raise ValueError(f"nEvents All={processed} is outside (0, {args.max_events}]")
        if counts[0] > processed:
            raise ValueError("Signal veto checks exceed processed events")
        # FilledHistograms bit 0 enables event-stage filling (the runner uses 31).
        if int(values.get("FilledHistograms", "0")) & 1:
            vz_stage = float(events.GetBinContent(9))
            if vz_stage != counts[0]:
                raise ValueError("nEvents vz stage differs from signal veto checks")
        selected = float(required("eventQA/eventActivity").GetEntries())
        jet_entries = float(required("inclusiveJet").GetEntries())
        eec_entries = float(required("energyEnergyCorrelator").GetEntries())
        report["statistics"] = {"processed_events": processed, "selected_events": selected,
                                "jet_entries": jet_entries, "eec_entries": eec_entries}
        if selected != counts[1]:
            raise ValueError("Selected event entries differ from signal veto kept count")
        empty = [name for name, value in (("selected_events", selected), ("jet_entries", jet_entries), ("eec_entries", eec_entries)) if value <= 0]
        if empty:
            report["empty_physics_statistics"] = empty
            if not args.allow_empty_physics:
                raise ValueError("Empty physics statistics: " + ", ".join(empty) + "; increase the event cap or use --allow-empty-physics")
            report["status"] = "INCONCLUSIVE"
        else:
            report["status"] = "PASS_BOUNDED_MC_SOFTWARE_CHECK"
    finally:
        source.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--card", type=Path, required=True)
    parser.add_argument("--max-events", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-empty-physics", action="store_true")
    args = parser.parse_args()
    if args.max_events <= 0:
        parser.error("--max-events must be positive")
    report = {"status": "FAIL", "root_file": str(args.root_file), "card": str(args.card),
              "scope": "Bounded MC software validation; no campaign coverage or physics closure certification"}
    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        check_file(ROOT, args, report)
    except Exception as error:
        report["error"] = str(error)
        report["status"] = "FAIL"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"{report['status']}: {args.output}")
    if "error" in report:
        print(report["error"])
    return int(report["status"] == "FAIL")


if __name__ == "__main__":
    raise SystemExit(main())
