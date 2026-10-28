#!/usr/bin/env python3
"""Reduce retained holdout event vectors into joint raw-yield statistics.

The covariance is sum_e x_e x_e^T (Poisson event resampling, fixed maps).
This is a diagnostic. It excludes map/response MC uncertainty, mixed-event
reuse, normalized spectra, unfolding and systematics; those remain NOT_RUN.
"""
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

COMPONENTS = ("gen", "raw", "single", "pair_supported", "masked_single")


def reduce_reports(reports):
    folds, events = {}, defaultdict(dict)
    for report in reports:
        for row in report["event_folds"]:
            key = tuple(row["id"])
            value = row["fold"]
            if key in folds and folds[key] != value:
                raise ValueError("Training/validation event leakage across reports")
            folds[key] = value
        for row in report["holdout_event_vectors"]:
            event = tuple(row["id"])
            if row["fold"] != 1 or folds.get(event) != 1:
                raise ValueError("Holdout vector does not have a held-out physical identity")
            events[event]  # Retain zero-contribution holdout events in exposure.
            for channel in row["channels"]:
                for component in COMPONENTS:
                    column = (*channel["key"], component)
                    if column in events[event]:
                        raise ValueError("Duplicate event/channel: overlapping inputs must not be counted twice")
                    value = float(channel[component])
                    if not math.isfinite(value):
                        raise ValueError("Nonfinite event sufficient statistic")
                    events[event][column] = value
    columns = sorted({column for vector in events.values() for column in vector})
    indices = {column: i for i, column in enumerate(columns)}
    total, covariance = [0.] * len(columns), defaultdict(float)
    for vector in events.values():
        values = sorted((indices[column], value) for column, value in vector.items() if value)
        for pos, (i, first) in enumerate(values):
            total[i] += first
            for j, second in values[pos:]:
                covariance[i, j] += first * second
    return dict(status="DIAGNOSTIC_FIXED_MAP_HOLDOUT_ONLY",
                covariance_definition="sum of physical-event outer products; Poisson(1) event resampling",
                columns=[list(column) for column in columns], totals=total,
                covariance_upper=[[i, j, value] for (i, j), value in sorted(covariance.items())],
                physical_events=len(events), contributing_physical_events=sum(bool(v) for v in events.values()), event_ids=[list(key) for key in sorted(events)],
                finite_mc_map_uncertainty="NOT_RUN", mixed_event_reuse="NOT_RUN",
                final_chain_covariance="NOT_RUN",
                pair_supported_warning="partial supported-map yield; masked_single is retained; not a completed corrected spectrum")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = reduce_reports([json.loads(path.read_text()) for path in args.reports])
    result["input_reports"] = [str(path) for path in args.reports]
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(result["status"], result["physical_events"], "physical holdout events")


if __name__ == "__main__":
    main()
