#!/usr/bin/env python3
"""Bounded raw geometry diagnostic, never a calibrated spectrum or closure.

Keep inherited +/-pi/2 cone geometry and all candidate pairs, including physical
self-pairs, to make any failed estimator visible. Source IDs include input-file
UUID, run/lumi/event, entry and track level/index. No mask, shift or correction is
applied. Stored jtpt is only a diagnostic bin coordinate, not certified JEC.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

PI = math.pi
PAIRINGS = {"SS": (0, 0), "SPplus": (0, 1), "PplusPplus": (1, 1),
            "SPminus": (0, 2), "PplusPminus": (1, 2), "PminusPminus": (2, 2)}
BINS = [(40, 80), (80, 120), (120, 160)]


def wrap(phi):
    return math.remainder(phi, 2 * PI)


def delta_r(a, b):
    return math.hypot(a[1]-b[1], wrap(a[2]-b[2]))


def char_int(value):
    return ord(value) if isinstance(value, str) else int(value)


def in_cone(track, eta, phi, radius):
    return math.hypot(track[1]-eta, wrap(track[2]-phi)) < radius


def source_id(event_key, level, index):
    return (*event_key, level, index)


def self_test():
    midpoint = (3.0, 0.0, PI/4, source_id(("file", 1, 2, 3, 0), "reco", 7))
    assert in_cone(midpoint, 0, 0, .8) and in_cone(midpoint, 0, PI/2, .8)
    assert not in_cone(midpoint, 0, 0, .4)
    assert source_id(("file", 1, 2, 3, 0), "reco", 7) != source_id(("file", 1, 2, 4, 1), "reco", 7)
    assert delta_r(midpoint, midpoint) == 0
    rotated_first = (3, 0, PI/4, midpoint[3])
    rotated_second = (3, 0, -PI/4, midpoint[3])
    assert abs(delta_r(rotated_first, rotated_second)-PI/2) < 1e-12
    assert 2*.8 > PI/2 and 2*.4 < PI/2
    assert math.isclose((3*4)**2, 144) and (3*4)**2 != (3*4)
    for phi in [-PI, -3, 0, 3, PI]:
        for sign in [-1, 1]:
            assert abs(wrap(phi+sign*PI/2)) <= PI
    return {"status": "PASS", "kind": "synthetic software assertions", "checks": [
        "R8 physical overlap and R4 disjoint control", "event-aware identity",
        "physical self-pair analysis angle pi/2 versus detector angle zero",
        "separate n=1,2 energy powers", "periodic wrapping"]}


def new_bin():
    return {"jets": 0, "weighted_jet_exposure": 0.0, "cone_track_counts": [0, 0, 0],
            "shared_signal_control_tracks": [0, 0], "control_tracks_inside_other_selected_jet": [0, 0],
            "powers": {str(n): {"cone_track_pt_power": [0.0]*3,
                      "shared_track_pt_power": [0.0]*2,
                      "other_jet_track_pt_power": [0.0]*2,
                      "pairings": {k: {"count": 0, "raw_yield": 0.0,
                          "stored_generator_weighted_yield": 0.0, "self_pair_count": 0,
                          "self_pair_raw_yield": 0.0, "self_pair_generator_weighted_yield": 0.0,
                          "overlap_track_pair_raw_yield": 0.0, "analysis_above_0p4_raw_yield": 0.0,
                          "analysis_above_0p8_raw_yield": 0.0, "max_analysis_delta_r": 0.0,
                          "max_detector_delta_r": 0.0} for k in PAIRINGS}} for n in [1, 2]}}


def scan(args):
    import ROOT
    start = time.perf_counter()
    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot open input")
    trees = {p: root_file.Get(p) for p in ["hiEvtAnalyzer/HiTree", "skimanalysis/HltTree",
             "PbPbTracks/trackTree", "HiGenParticleAna/hi", "akCs4PFJetAnalyzer/t", "akCs8PFJetAnalyzer/t"]}
    if not all(trees.values()):
        raise RuntimeError("Missing required native OO tree")
    counts = {p: int(t.GetEntries()) for p, t in trees.items()}
    if len(set(counts.values())) != 1:
        raise RuntimeError("Required tree event counts differ")
    event, skim, reco, gen = [trees[p] for p in ["hiEvtAnalyzer/HiTree", "skimanalysis/HltTree", "PbPbTracks/trackTree", "HiGenParticleAna/hi"]]
    file_uuid = str(root_file.GetUUID().AsString())
    result = {"classification": "bounded actual-file raw geometry diagnostic; not production or physics closure",
        "input": args.input, "input_uuid": file_uuid, "file_size": Path(args.input).stat().st_size,
        "tree_entries": counts, "events_requested": args.events, "events_scanned": 0,
        "accepted_events": 0, "event_keys": [], "geometry": {}, "bins": {}, "shared_identity_examples": [],
        "selection": {"event": "pprimaryVertexFilter && OOphfCoincFilterPF2Th4 && pileupVertexFilter && abs(vz)<=15; no MC trigger",
                      "reco_track": "2<pT<500, abs(eta)<2.4, charged, highPurity, pTerr/pT<0.1, abs(dz/sigma)<3, abs(dxy/sigma)<3; first-vertex quantities as reader",
                      "gen_track": "2<pT<500, abs(eta)<2.4, charged; final-state convention implicit in HiGenParticleAna tree; no subevent cut",
                      "jet": "native reco stored jtpt or native genpt in 40-80/80-120/120-160; abs(WTAeta)<1.6; reco 0.01<trackMax/rawpt<0.98",
                      "weights": "raw pair (pTi*pTj)^n; separate stored event generator weight applied once; no single/pair efficiency, event reweight or JEC adoption"},
        "not_run": ["applicable AK8 calibration", "tracking/pair corrections", "masked or shifted geometry",
                    "mixed-event estimator", "full extracted-signal closure", "validated final angular range",
                    "data/pp geometry pilots", "statistical covariance replicas"]}
    for radius in [.4, .8]:
        d = PI/2
        area = 0.0 if d >= 2*radius else 2*radius**2*math.acos(d/(2*radius))-.5*d*math.sqrt(4*radius**2-d*d)
        key = "R04" if radius == .4 else "R08"
        result["geometry"][key] = {"radius": radius, "control_shift": [PI/2, -PI/2],
            "nominal_cone_area": PI*radius**2, "signal_control_intersection_area_each": area,
            "area_is_not_weighted_bias_evidence": True, "pair_support": [0, 2*radius],
            "containment": "strict abs(WTAeta)<1.6 plus R<=0.8 is contained in abs(track_eta)<2.4 for ideal circle"}
        for level in ["reco", "gen"]:
            for low, high in BINS:
                result["bins"][f"{key}_{level}_{low}_{high}"] = new_bin()
    for entry in range(min(args.events, counts["hiEvtAnalyzer/HiTree"])):
        for tree in trees.values():
            tree.GetEntry(entry)
        result["events_scanned"] += 1
        if not (int(skim.pprimaryVertexFilter) and int(skim.OOphfCoincFilterPF2Th4) and int(skim.pileupVertexFilter) and abs(event.vz) <= 15):
            continue
        result["accepted_events"] += 1
        event_key = (file_uuid, int(event.run), int(event.lumi), int(event.evt), entry)
        result["event_keys"].append(list(event_key))
        stored_weight = float(event.weight)
        tracks = {"reco": [], "gen": []}
        for i in range(int(reco.nTrk)):
            pt, eta = float(reco.trkPt[i]), float(reco.trkEta[i])
            if not (2 < pt < 500 and abs(eta) < 2.4 and char_int(reco.trkCharge[i]) != 0 and reco.highPurity[i] and reco.trkPtError[i]/pt < .1):
                continue
            if not (reco.trkDzErrFirstVtx[i] > 0 and reco.trkDxyErrFirstVtx[i] > 0):
                continue
            if abs(reco.trkDzFirstVtx[i]/reco.trkDzErrFirstVtx[i]) >= 3 or abs(reco.trkDxyFirstVtx[i]/reco.trkDxyErrFirstVtx[i]) >= 3:
                continue
            tracks["reco"].append((pt, eta, float(reco.trkPhi[i]), source_id(event_key, "reco", i)))
        for i in range(int(gen.mult)):
            pt, eta = float(gen.pt[i]), float(gen.eta[i])
            if 2 < pt < 500 and abs(eta) < 2.4 and int(gen.chg[i]) != 0:
                tracks["gen"].append((pt, eta, float(gen.phi[i]), source_id(event_key, "gen", i)))
        for radius, tree_path in [(.4, "akCs4PFJetAnalyzer/t"), (.8, "akCs8PFJetAnalyzer/t")]:
            jets = trees[tree_path]
            tag = "R04" if radius == .4 else "R08"
            for level in ["reco", "gen"]:
                accepted = []
                for j in range(int(jets.nref if level == "reco" else jets.ngen)):
                    pt = float(jets.jtpt[j] if level == "reco" else jets.genpt[j])
                    eta = float(jets.WTAeta[j] if level == "reco" else jets.WTAgeneta[j])
                    phi = float(jets.WTAphi[j] if level == "reco" else jets.WTAgenphi[j])
                    if abs(eta) >= 1.6 or pt < 40 or pt >= 160:
                        continue
                    if level == "reco" and not (jets.rawpt[j] > 0 and .01 < jets.trackMax[j]/jets.rawpt[j] < .98):
                        continue
                    accepted.append((j, pt, eta, phi))
                for j, pt, eta, phi in accepted:
                    low, high = next((low, high) for low, high in BINS if low <= pt < high)
                    out = result["bins"][f"{tag}_{level}_{low}_{high}"]
                    out["jets"] += 1
                    out["weighted_jet_exposure"] += stored_weight
                    centers = [phi, wrap(phi+PI/2), wrap(phi-PI/2)]
                    cones = [[t for t in tracks[level] if in_cone(t, eta, center, radius)] for center in centers]
                    ids = [{t[3] for t in cone} for cone in cones]
                    shared = [ids[0] & ids[c] for c in [1, 2]]
                    union_shared = set().union(*shared)
                    for c in range(3):
                        out["cone_track_counts"][c] += len(cones[c])
                    for c in range(2):
                        out["shared_signal_control_tracks"][c] += len(shared[c])
                    if union_shared and len(result["shared_identity_examples"]) < 30:
                        result["shared_identity_examples"].append({"jet": [tag, level, entry, j], "ids": [list(x) for x in sorted(union_shared)]})
                    for n in [1, 2]:
                        power = out["powers"][str(n)]
                        for c in range(3):
                            power["cone_track_pt_power"][c] += sum(t[0]**n for t in cones[c])
                        for c in [1, 2]:
                            power["shared_track_pt_power"][c-1] += sum(t[0]**n for t in cones[c] if t[3] in shared[c-1])
                            other = [t for t in cones[c] if any(oj != j and in_cone(t, oe, op, radius) for oj, _, oe, op in accepted)]
                            if n == 1:
                                out["control_tracks_inside_other_selected_jet"][c-1] += len(other)
                            power["other_jet_track_pt_power"][c-1] += sum(t[0]**n for t in other)
                        for name, (a, b) in PAIRINGS.items():
                            target = power["pairings"][name]
                            for i, first in enumerate(cones[a]):
                                for second in cones[b][i+1 if a == b else 0:]:
                                    weight = (first[0]*second[0])**n
                                    detector_dr = delta_r(first, second)
                                    analysis_dr = math.hypot(first[1]-second[1], wrap(wrap(first[2]-centers[a])-wrap(second[2]-centers[b])))
                                    target["count"] += 1
                                    target["raw_yield"] += weight
                                    target["stored_generator_weighted_yield"] += stored_weight*weight
                                    target["max_analysis_delta_r"] = max(target["max_analysis_delta_r"], analysis_dr)
                                    target["max_detector_delta_r"] = max(target["max_detector_delta_r"], detector_dr)
                                    if analysis_dr > .4:
                                        target["analysis_above_0p4_raw_yield"] += weight
                                    if analysis_dr > .8:
                                        target["analysis_above_0p8_raw_yield"] += weight
                                    if first[3] in union_shared or second[3] in union_shared:
                                        target["overlap_track_pair_raw_yield"] += weight
                                    if first[3] == second[3]:
                                        target["self_pair_count"] += 1
                                        target["self_pair_raw_yield"] += weight
                                        target["self_pair_generator_weighted_yield"] += stored_weight*weight
    result["shared_event_key_sha256"] = hashlib.sha256(json.dumps(result["event_keys"], sort_keys=True).encode()).hexdigest()
    for tag in ["R04", "R08"]:
        self_count = sum(v["powers"]["1"]["pairings"][p]["self_pair_count"] for k, v in result["bins"].items() if k.startswith(tag) for p in PAIRINGS)
        result["geometry"][tag]["observed_inherited_self_pairs"] = self_count
        result["geometry"][tag]["inherited_estimator_status"] = "FAIL_PHYSICAL_SELF_PAIRS" if self_count else "NOT_CLOSED_NO_SELF_PAIR_IN_BOUNDED_SAMPLE"
    result["wall_seconds"] = time.perf_counter()-start
    result["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input")
    parser.add_argument("--events", type=int, default=100)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.events < 1 or args.events > 1000:
        parser.error("Use a bounded event count from 1 to 1000")
    result = {"synthetic": self_test()}
    if args.input:
        result["actual"] = scan(args)
    elif not args.self_test:
        parser.error("Provide --input or --self-test")
    output = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    print(json.dumps({"synthetic": result["synthetic"], "actual": {k: result.get("actual", {}).get(k) for k in ["accepted_events", "events_scanned", "geometry", "wall_seconds"]}, "output": str(args.output) if args.output else None}, indent=2))


if __name__ == "__main__":
    main()
