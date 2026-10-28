#!/usr/bin/env python3
"""Primary channel metadata gate. No physics inputs or validation are supplied."""
import argparse
import json
import math

JET_BINS = ((40, 80), (80, 120), (120, 160))
RADII = (0.4, 0.8)
POWERS = (1, 2)


def channels():
    return [dict(channel=f"INC_R{round(r * 10):02d}_J{lo}_{hi}_N{n}",
                 radius=r, jet_pt=[lo, hi], energy_power=n,
                 production_status="NOT_RUN", physics_status="NOT_RUN")
            for r in RADII for lo, hi in JET_BINS for n in POWERS]


def validate_ratio(oo, pp):
    """Reject incompatible metadata before dividing primary spectra."""
    common = ("energy_tev", "radius", "axis", "jet_pt", "track_pt_min",
              "track_eta_max", "track_cut_version", "angular_edges",
              "normalization_range", "stage", "energy_power")
    for field in common:
        if field not in oo or field not in pp or oo[field] is None or pp[field] is None:
            raise ValueError(f"Missing validated ratio metadata: {field}")
        if oo[field] != pp[field]:
            raise ValueError(f"Incompatible ratio metadata: {field}")
    if oo.get("system") != "OO" or pp.get("system") != "pp":
        raise ValueError("Ratio requires OO numerator and pp denominator")
    if oo["radius"] not in RADII or tuple(oo["jet_pt"]) not in JET_BINS:
        raise ValueError("Channel is outside the primary measurement contract")
    if oo["energy_power"] not in POWERS or oo["energy_tev"] != 5.36 or oo["axis"] != "WTA":
        raise ValueError("Unsupported energy, axis or energy power")
    if oo["track_pt_min"] != 2 or oo["track_eta_max"] != 2.4:
        raise ValueError("Primary track cuts differ from the contract")
    edges, limits = oo["angular_edges"], oo["normalization_range"]
    if (not isinstance(edges, list) or len(edges) < 2 or
            not all(isinstance(x, (float, int)) and math.isfinite(x) for x in edges) or
            any(a >= b for a, b in zip(edges, edges[1:]))):
        raise ValueError("Invalid angular edges")
    if (not isinstance(limits, list) or len(limits) != 2 or
            not all(isinstance(x, (float, int)) and math.isfinite(x) for x in limits) or
            not edges[0] <= limits[0] < limits[1] <= edges[-1]):
        raise ValueError("Invalid normalization range")
    for sample in (oo, pp):
        for field in ("native_collection_verified", "correction_applicability_verified",
                      "independent_closure_verified", "covariance_verified",
                      "normalization_range_verified"):
            if sample.get(field) is not True:
                raise ValueError(f"BLOCKED: {sample['system']} {field}")
        if sample.get("stage") != "final_normalized":
            raise ValueError("Primary ratio requires the documented final correction chain")
        if not sample.get("input_hashes") or not sample.get("calibration_hashes"):
            raise ValueError("Missing input/calibration provenance")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-ratio", nargs=2, metavar=("OO_METADATA", "PP_METADATA"))
    args = parser.parse_args()
    if args.check_ratio:
        with open(args.check_ratio[0]) as stream:
            oo = json.load(stream)
        with open(args.check_ratio[1]) as stream:
            pp = json.load(stream)
        try:
            validate_ratio(oo, pp)
        except ValueError as error:
            parser.exit(2, str(error) + "\n")
        print("PASS: metadata gate only; inspect the referenced validation evidence")
    else:
        print(json.dumps(channels(), indent=2))


if __name__ == "__main__":
    main()
