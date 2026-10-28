#!/usr/bin/env python3
"""Reproducible bounded workflow software checks and optional real-file pilots."""
import argparse
import hashlib
import json
import shlex
import subprocess
import sys
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--oo", type=Path)
    parser.add_argument("--pp", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    checks = []

    def run(name, command, timeout=180):
        started = time.monotonic()
        with (out/(name+".log")).open("w") as log:
            log.write(shlex.join(command)+"\n")
            log.flush()
            result = subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        checks.append(dict(name=name, command=command, exit_code=result.returncode,
                           elapsed_seconds=time.monotonic()-started))
        (out/"commands.json").write_text(json.dumps(checks, indent=2)+"\n")
        print(name, "PASS" if result.returncode == 0 else "FAIL", flush=True)
        if result.returncode:
            raise RuntimeError(name+" failed; see "+str(out/(name+".log")))

    run("build", ["make", "-B", "-j4"])
    run("cards", [sys.executable, "validation/run3Workflow/check_cards.py"])
    run("metadata_statistics", [sys.executable, "-m", "unittest", "discover", "-s", "validation/run3Workflow", "-p", "test_*.py", "-v"])
    run("projection_synthetic", [sys.executable, "validation/run3Workflow/check_projection.py", str(out/"projection-fixtures")])
    flags = shlex.split(subprocess.check_output(["root-config", "--cflags", "--libs"], text=True))
    objects = sorted(str(p.relative_to(repo)) for p in (repo/"src").glob("*.o"))
    binary = out/"checkAnalyzerConfiguration"
    run("build_analyzer_regression", ["g++", "-fno-access-control", "-I.",
         "validation/run3Response/checkAnalyzerConfiguration.cxx", *objects, *flags, "-lEG", "-o", str(binary)])
    run("analyzer_regression", [str(binary), str(out/"analyzer-fixtures")])
    for system, path in (("OO", args.oo), ("pp", args.pp)):
        if path is None:
            continue
        target = out/(system+"_R04_raw.root")
        card = "validation/run3Workflow/cards/"+system+"MC_R04_diagnostic.input"
        run(system+"_actual_capped", [str(repo/"eecAnalysis"), str(path.resolve()), card, str(target), "1", "true"], timeout=120)
        projected = out/(system+"_R04_projected.root")
        invocation = "validation/run3Workflow/projectInclusiveRaw.C("+json.dumps(str(target))+","+json.dumps(str(projected))+")"
        run(system+"_actual_projection", ["root", "-l", "-b", "-q", invocation])
    summary = dict(status="PASS_SOFTWARE_AND_BOUNDED_DIAGNOSTICS_ONLY", physics_closure="NOT_RUN",
                   commands=checks, head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
                   outputs={str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in out.rglob("*") if p.is_file() and p.suffix in (".root", ".log", ".json")})
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
