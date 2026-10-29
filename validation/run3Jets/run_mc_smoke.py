#!/usr/bin/env python3
"""Run bounded OO/ppRef MC checks using scratch copies of production cards.

Run from a configured CMSSW/ROOT environment. Existing output directories are
refused. No production card, input forest, or batch submission is modified.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys


REPO = Path(__file__).resolve().parents[2]
DEFAULT_INPUTS = {
    "OO": "/eos/user/j/jvelazqu/Foresting/OOMC/unembedded_L1HLT_objects_MC_1.root",
    "pp": "/eos/user/j/jvelazqu/Foresting/ppMC/HiForestMiniAOD_MC_1.root",
}
BASE_CARDS = {"OO": "cardEECOOMC.input", "pp": "cardEECpp5v36TeVMC.input"}


def card_values(text):
    """Parse active keys and fail on duplicates rather than shadowing them."""
    values = {}
    for number, line in enumerate(text.splitlines(), 1):
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        if len(fields) < 2:
            raise ValueError(f"Card line {number} has no value: {line}")
        if fields[0] in values:
            raise ValueError(f"Duplicate card key: {fields[0]}")
        values[fields[0]] = " ".join(fields[1:])
    return values


def bounded_card(text, events, correlation_type):
    values = card_values(text)
    for key in ("Run3JetID", "ApplyJetVeto"):
        if values.get(key) != "1":
            raise ValueError(f"Production card must already set {key} 1")
    for key in ("JetVetoMapFile", "JetVetoMapHistogram", "JetVetoPtConvention", "JetVetoInputMinPt"):
        if key not in values:
            raise ValueError(f"Production card is missing {key}")
    overrides = {
        "MaxEventsPerFile": str(events),
        "BackgroundMethods": "0",
        "McCorrelationType": str(correlation_type),
        "TriggerSelection": "0",
        "FilledHistograms": "31",
        "SkipCovarianceMatrix": "1",
        "DebugLevel": "1",
    }
    lines = ["# Bounded MC software test; overrides are local to this scratch card."]
    for line in text.splitlines():
        fields = line.split("#", 1)[0].split()
        if fields and fields[0] in overrides:
            key = fields[0]
            lines.append(f"{key} {overrides.pop(key)} # bounded MC smoke test")
        else:
            lines.append(line)
    lines.extend(f"{key} {value} # bounded MC smoke test" for key, value in overrides.items())
    result = "\n".join(lines) + "\n"
    scratch_values = card_values(result)
    for key in ("Run3JetID", "ApplyJetVeto", "JetVetoMapFile", "JetVetoMapHistogram", "JetVetoPtConvention", "JetVetoInputMinPt"):
        if scratch_values[key] != values[key]:
            raise ValueError(f"Smoke preparation unexpectedly changed {key}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=("OO", "pp", "both"), default="both")
    parser.add_argument("--events", type=int, default=1000, help="Maximum events from the one input file per system")
    parser.add_argument("--output", type=Path, required=True, help="New scratch directory; must not already exist")
    parser.add_argument("--oo-input", default=DEFAULT_INPUTS["OO"])
    parser.add_argument("--pp-input", default=DEFAULT_INPUTS["pp"])
    parser.add_argument("--correlation-type", type=int, choices=(0, 3), default=0, help="0=RecoReco (default), 3=GenGen with reconstructed event veto")
    parser.add_argument("--prepare-only", action="store_true", help="Write scratch cards, manifest, and commands without running the analyzer")
    parser.add_argument("--allow-empty-physics", action="store_true", help="Report empty jet/EEC statistics as INCONCLUSIVE instead of failing")
    args = parser.parse_args()
    if args.events <= 0:
        parser.error("--events must be positive; uncapped tests are not supported")
    output = args.output.expanduser().resolve()
    if output.exists():
        parser.error(f"Output already exists; choose a fresh directory: {output}")
    analyzer = REPO / "eecAnalysis"
    if not analyzer.is_file():
        parser.error("Build ./eecAnalysis first (make -j4)")
    systems = ("OO", "pp") if args.system == "both" else (args.system,)
    prepared = []
    for system in systems:
        input_file = getattr(args, "oo_input" if system == "OO" else "pp_input")
        if not input_file.startswith("root://") and not Path(input_file).is_file():
            parser.error(f"Input file does not exist: {input_file}")
        base = REPO / BASE_CARDS[system]
        original = base.read_text()
        prepared.append((system, input_file, base, original, bounded_card(original, args.events, args.correlation_type)))
    output.mkdir(parents=True)
    jobs = []
    commands = ["#!/usr/bin/env bash", "set -euo pipefail", f"cd {shlex.quote(str(REPO))}"]
    for system, input_file, base, original, scratch in prepared:
        stem = f"{system}_R04_{'RecoReco' if args.correlation_type == 0 else 'GenGen'}"
        card = output / f"{stem}.input"
        card.write_text(scratch)
        root_file = output / f"{stem}.root"
        log = output / f"{stem}.log"
        command = [str(analyzer), input_file, str(card), str(root_file), "1", "true"]
        check = [sys.executable, str(REPO / "validation/run3Jets/check_mc_smoke.py"), str(root_file), "--card", str(card), "--max-events", str(args.events), "--output", str(output / f"{stem}_check.json")]
        if args.allow_empty_physics:
            check.append("--allow-empty-physics")
        jobs.append({"system": system, "input": input_file, "base_card": str(base), "base_card_sha256": hashlib.sha256(original.encode()).hexdigest(), "scratch_card": str(card), "output": str(root_file), "log": str(log), "command": command, "check_command": check, "status": "PREPARED"})
        commands += [f"{shlex.join(command)} > {shlex.quote(str(log))} 2>&1", shlex.join(check)]
    manifest = {"purpose": "BOUNDED_MC_SOFTWARE_CHECK; no production coverage/physics closure certification", "max_events_per_file": args.events, "jobs": jobs}
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "commands.sh").write_text("\n".join(commands) + "\n")
    print(f"Scratch cards and reproducible commands: {output}", flush=True)
    if args.prepare_only:
        print(f"Run: bash {shlex.quote(str(output / 'commands.sh'))}")
        return
    for job in jobs:
        print(f"Running {job['system']} (maximum {args.events} events); log: {job['log']}", flush=True)
        with Path(job["log"]).open("w") as log:
            result = subprocess.run(job["command"], cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
        job["analyzer_exit_code"] = result.returncode
        job["status"] = "ANALYZER_FAILED" if result.returncode else "ANALYZER_COMPLETE"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        if result.returncode:
            raise SystemExit(f"Analyzer failed; inspect {job['log']}")
        checked = subprocess.run(job["check_command"], cwd=REPO)
        job["check_exit_code"] = checked.returncode
        job["status"] = "CHECK_FAILED" if checked.returncode else "CHECK_COMPLETE"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        if checked.returncode:
            raise SystemExit(f"Output check failed for {job['system']}")
    print(f"Finished bounded MC checks. Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
