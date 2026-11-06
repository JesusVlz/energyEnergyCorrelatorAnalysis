#!/usr/bin/env python3
"""Run a bounded native OO closure pilot and save exact executable provenance.

Build eecAnalysis first. This runner never expands caps, submits jobs, or
reuses a ROOT/ledger output prefix. The source card is copied with explicit
per-run global caps and an output prefix; production cards are unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import tarfile
import time

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CARD = REPO / "validation/run3Mixing/cards/OO_GenGen_R04_closure.input"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def capture(command):
    result = subprocess.run(command, cwd=REPO, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=True)
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=REPO / "eecAnalysis")
    parser.add_argument("--signal", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--card", type=Path, default=DEFAULT_CARD)
    parser.add_argument("--signal-events", type=int, default=10000)
    parser.add_argument("--pool-events", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20261030)
    args = parser.parse_args()
    if not 0 < args.signal_events <= 10000 or not 0 < args.pool_events <= 100000:
        parser.error("This bounded pilot permits at most 10000 signal / 100000 total pool events")
    if args.seed <= 0 or args.seed > 2**31-1:
        parser.error("Seed must be a positive exactly representable card integer")
    if int(float(args.seed)) != args.seed:
        parser.error("Invalid seed")
    import struct
    if int(struct.unpack("f", struct.pack("f", args.seed))[0]) != args.seed:
        parser.error("Seed is not exactly representable by the card Float_t")
    signal = args.signal.resolve(strict=True)
    executable = args.executable.resolve(strict=True)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    prefix = out / "pilot"
    protected = [out / name for name in ("pilot.root", "pilot.input", "run.json", "run.log", "commands.sh", "source.diff", "source_snapshot.tar.gz", "signal.txt")]
    protected += list(out.glob("pilot.*.csv"))
    if any(path.exists() for path in protected):
        raise FileExistsError("Choose a fresh run prefix: " + str(prefix))
    text = args.card.read_text()
    overrides = dict(OOClosureMaxSignalEvents=args.signal_events,
                     OOClosureMaxPoolEvents=args.pool_events,
                     OOClosureSeed=args.seed, OOClosureOutputPrefix=str(prefix))
    for key, value in overrides.items():
        text, count = re.subn(r"^"+key+r"\s+.*$", key+" "+str(value), text, flags=re.M)
        if count != 1:
            raise ValueError("Expected exactly one " + key)
    card = out / "pilot.input"
    card.write_text(text)
    (out / "signal.txt").write_text(str(signal)+"\n")
    command = [str(executable), str(signal), str(card), str(out / "pilot.root"), "1", "true"]
    (out / "commands.sh").write_text("#!/bin/bash\nset -euo pipefail\ncd "+shlex.quote(str(REPO))+"\n"+shlex.join(command)+"\n")
    (out / "source.diff").write_text(capture(["git", "diff", "--binary", "HEAD"]))
    sources = [REPO / x for x in capture(["git", "ls-files", "--cached", "--others", "--exclude-standard"]).splitlines()
               if x.startswith(("src/", "plotting/", "validation/run3Mixing/"))
               and "/results/" not in x and Path(x).suffix in (".h", ".cxx", ".C", ".py", ".input")]
    sources += [REPO / "Makefile", REPO / "eecAnalysis"]
    with tarfile.open(out / "source_snapshot.tar.gz", "x:gz") as archive:
        for path in sources:
            if path.is_file() and path.name != "eecAnalysis":
                archive.add(path, arcname=str(path.relative_to(REPO)))
    pool_line = next(line.split("#")[0].split(maxsplit=1)[1].strip()
                     for line in text.splitlines() if line.startswith("MixingFileList "))
    pool_list = REPO / pool_line
    input_manifest = pool_list.with_suffix(".manifest.json")
    manifest = dict(status="RUNNING", started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    command=command, executable=str(executable), executable_sha256=sha(executable), source_commit=capture(["git", "rev-parse", "HEAD"]),
                    source_status=capture(["git", "status", "--short"]),
                    source_sha256={str(p.relative_to(REPO)): sha(p) for p in sources if p.is_file()},
                    source_snapshot_sha256=sha(out / "source_snapshot.tar.gz"),
                    root_version=capture(["root-config", "--version"]),
                    compiler_version=capture(["g++", "--version"]).splitlines()[0],
                    signal_file=str(signal), signal_file_bytes=signal.stat().st_size,
                    pool_list=str(pool_list), pool_list_sha256=sha(pool_list),
                    pool_schema_manifest=str(input_manifest),
                    pool_schema_manifest_sha256=sha(input_manifest) if input_manifest.is_file() else None,
                    card_sha256=sha(card), global_signal_cap=args.signal_events,
                    global_pool_cap=args.pool_events, seed=args.seed,
                    matching=dict(nh="stable charged HIJING pT>0.7 GeV, abs(eta)<2.4", dz_max_cm=.5,
                                  dnh_max="max(1,0.01*NHsignal)", hf_transform="OFF; exact raw identity"),
                    calibration_sha256={str(p.relative_to(REPO)): sha(p) for p in (
                        REPO/"jetEnergyCorrections/L2Relative_AK4PF_OO_V1.txt",
                        REPO/"trackCorrectionTables/OO2025/Eff_OO_2025_Hijing_MB_Centrality_fromHihfpf_NoPU_3D_Nominal_Official_18Nov2025.root")},
                    corrections="GenGen physical pT; reconstruction single-track/pair unity; Run-3 payloads retained for reco paths",
                    background_independence="Observed raw identities and nonempty stored-particle fingerprints globally disjoint; generator lineage unavailable",
                    uncertainty="Conditional fixed-pool joint physical-signal-event moments; donor reuse is not independent exposure")
    record = out / "run.json"
    record.write_text(json.dumps(manifest, indent=2)+"\n")
    started = time.monotonic()
    with (out / "run.log").open("x") as log:
        log.write(shlex.join(command)+"\n")
        log.flush()
        result = subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
    manifest.update(status="PASS_NATIVE_EXECUTION" if result.returncode == 0 else "FAILED", exit_code=result.returncode,
                    elapsed_seconds=time.monotonic()-started,
                    completed_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    manifest["output_sha256"] = {p.name: sha(p) for p in out.glob("pilot.*") if p.is_file()}
    record.write_text(json.dumps(manifest, indent=2)+"\n")
    print(manifest["status"], "seconds=", round(manifest["elapsed_seconds"], 2), "log=", out / "run.log", flush=True)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
