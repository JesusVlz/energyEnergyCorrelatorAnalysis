#!/usr/bin/env python3
"""Run synthetic checks and a bounded real OO diagnostic with separate logs."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oo", required=True)
    parser.add_argument("--pp", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events", type=int, default=100)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    folder = Path(__file__).resolve().parent
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    commands = []
    def run(name, command):
        result = subprocess.run(command, cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output/(name+".log")).write_text(shlex.join(command)+"\n"+result.stdout)
        commands.append({"name": name, "command": command, "exit_code": result.returncode})
        print(name+": "+("PASS" if result.returncode == 0 else "FAIL"), flush=True)
        if result.returncode:
            print(result.stdout)
            raise RuntimeError(name+" failed")
    run("buildGeometry", ["g++", "-std=c++17", "-Wall", "-Wextra", str(folder/"testRun3JetGeometry.cxx"), "-o", str(output/"testRun3JetGeometry")])
    run("geometrySynthetic", [str(output/"testRun3JetGeometry")])
    flags = shlex.split(subprocess.check_output(["root-config", "--cflags", "--libs"], text=True))
    run("buildBranchGuard", ["g++", "-std=c++17", "-Wall", str(folder/"testRun3JetBranches.cxx"), *flags, "-o", str(output/"testRun3JetBranches")])
    run("missingBranchSynthetic", [str(output/"testRun3JetBranches")])
    run("buildNativeReader", ["g++", "-std=c++17", "-Wall", str(folder/"testNativeGeneratorReader.cxx"),
        "src/ForestReader.cxx", "src/GeneratorLevelForestReader.cxx", *flags, "-o", str(output/"testNativeGeneratorReader")])
    run("nativeReaderActual", [str(output/"testNativeGeneratorReader"), args.oo, args.pp])
    run("geometryActual", [sys.executable, str(folder/"auditConeGeometry.py"), "--self-test", "--input", args.oo,
        "--events", str(args.events), "--output", str(output/f"coneGeometry{args.events}.json")])
    summary = {"status": "PASS_SOFTWARE_AND_BOUNDED_DIAGNOSTIC_ONLY", "physics_closure": "NOT_RUN",
               "commands": commands, "input_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
               "code_sha256": {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in
                   [folder/"testRun3JetGeometry.cxx", folder/"testRun3JetBranches.cxx", folder/"testNativeGeneratorReader.cxx", folder/"auditConeGeometry.py",
                    repo/"src/Run3JetGeometry.h", repo/"src/Run3JetBranchValidation.h", repo/"src/GeneratorLevelForestReader.cxx",
                    repo/"src/HighForestReader.cxx", repo/"src/UnfoldingForestReader.cxx", repo/"src/EECAnalyzer.cxx", repo/"Makefile", folder/"runChecks.py"]}}
    (output/"geometryChecks.json").write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
