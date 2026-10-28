#!/usr/bin/env python3
"""Prepare an immutable O-O or pp MC HTCondor campaign in AFS; never submit jobs."""

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile

try:
    from zoneinfo import ZoneInfo
except ImportError:
    ZoneInfo = None

MODES = {"RecoReco": 0, "RecoGen": 1, "GenReco": 2, "GenGen": 3, "Response": 3}
SAMPLES = {
    "OOMC": {"pattern": r"unembedded_L1HLT_objects_MC_([0-9]+)\.root",
             "card": "cardEECOOMC.input", "data_type": 10},
    "ppMC": {"pattern": r"HiForestMiniAOD_MC_([0-9]+)\.root",
             "card": "cardEECpp5v36TeVMC.input", "data_type": 12},
}
ARCH = "el9_amd64_gcc12"
CMSSW = "/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_0_11"


def timestamp():
    try:
        now = dt.datetime.now(ZoneInfo("America/Chicago")) if ZoneInfo else dt.datetime.now(dt.timezone.utc)
    except Exception:
        now = dt.datetime.now(dt.timezone.utc)
    return now


def absolute_path(value):
    path = Path(os.path.abspath(os.path.expanduser(str(value))))
    if any(char.isspace() for char in str(path)):
        raise ValueError("Paths must not contain whitespace: {}".format(path))
    if any(char in str(path) for char in '\"\'\n\r,;'):
        raise ValueError("Path contains unsupported submission syntax: {}".format(path))
    return path


def eos_url(path):
    if not str(path).startswith("/eos/"):
        raise ValueError("Expected an EOS path: {}".format(path))
    return "root://eosuser.cern.ch//" + str(path).lstrip("/")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_output(repo, *args):
    result = subprocess.run(["git", "-C", str(repo), *args], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return result.stdout.rstrip("\n")


def discover_inputs(input_dir, sample="OOMC"):
    input_pattern = re.compile(SAMPLES[sample]["pattern"])
    inputs = []
    for path in input_dir.iterdir():
        match = input_pattern.fullmatch(path.name)
        if match is None or not path.is_file():
            continue
        info = path.stat()
        if info.st_size == 0:
            raise ValueError("Input file is empty: {}".format(path))
        inputs.append({"index": int(match.group(1)), "path": str(path),
                       "url": eos_url(path), "stem": path.stem,
                       "size_bytes": info.st_size, "mtime_ns": info.st_mtime_ns,
                       "mtime_utc": dt.datetime.fromtimestamp(info.st_mtime, dt.timezone.utc).isoformat()})
    inputs.sort(key=lambda item: (item["index"], item["path"]))
    if not inputs:
        raise ValueError("No nonempty {} input forests in {}".format(sample, input_dir))
    return inputs


def card_text(base, mode, background_methods, sample="OOMC", match_eec_binning=False):
    overrides = {"DataType": str(SAMPLES[sample]["data_type"]), "McCorrelationType": str(MODES[mode]),
                 "TriggerSelection": "0", "MatchJets": "1" if mode == "Response" else "0",
                 "FilledHistograms": "479" if mode == "Response" else "31",
                 "JetType": "1", "JetAxis": "1", "JetRadius": "0.4",
                 "JetUncertainty": "0", "JetPtWeight": "0", "SmearDeltaR": "0",
                 "SmearEnergyWeight": "0", "DisableTrackPairEfficiencyCorrection": "1",
                 "SkipCovarianceMatrix": "1", "DebugLevel": "1",
                 "BackgroundMethods": "0" if mode == "Response" else str(background_methods)}
    if match_eec_binning:
        settings = {}
        for line in base.splitlines():
            fields = line.split("#", 1)[0].split()
            if fields:
                if fields[0] in settings:
                    raise ValueError("Duplicate card setting: {}".format(fields[0]))
                settings[fields[0]] = fields[1:]
        for key in ("JetPtBinEdgesEEC", "TrackPtBinEdgesEEC"):
            edges = [float(value) for value in settings.get(key, [])]
            if (len(edges) < 2 or not all(math.isfinite(edge) for edge in edges)
                    or any(not (left < right) for left, right in zip(edges, edges[1:]))):
                raise ValueError("{} must contain increasing bin edges".format(key))
        overrides["TrackPtBinEdges"] = " ".join(settings["TrackPtBinEdgesEEC"])
    found = set()
    lines = []
    for line in base.splitlines():
        stripped = line.lstrip()
        key = stripped.split(None, 1)[0] if stripped and not stripped.startswith("#") else None
        if key in overrides:
            if key in found:
                raise ValueError("Duplicate card setting: {}".format(key))
            found.add(key)
            comment = " #" + line.split("#", 1)[1] if "#" in line else ""
            line = "{} {}{}".format(key, overrides[key], comment)
        lines.append(line)
    missing = sorted(set(overrides) - found)
    if missing:
        raise ValueError("Base card is missing settings: {}".format(", ".join(missing)))
    if match_eec_binning:
        if "MatchEECBinning" in settings:
            raise ValueError("MatchEECBinning is controlled by --match-eec-binning; remove it from the base card")
        lines.append("MatchEECBinning 1 # Use EEC jet and track pT edges for sparse pT axes")
    return "# Prepared {} configuration: {}\n{}\n".format(sample, mode, "\n".join(lines))


def ignore_build_products(directory, names):
    return [name for name in names if name == "__pycache__" or name == ".git"
            or name.endswith((".o", ".a", ".so", ".d", ".pyc", ".dylib"))
            or name.endswith(".dSYM")]


def write_executable(path, content):
    path.write_text(content)
    path.chmod(0o755)


def build_snapshot(snapshot, cmssw_base, arch, jobs, log_path):
    commands = ["set -eo pipefail", "source /cvmfs/cms.cern.ch/cmsset_default.sh",
                "export SCRAM_ARCH={}".format(shlex.quote(arch)),
                "cd {}".format(shlex.quote(str(cmssw_base))),
                'eec_runtime=$(scram runtime -sh)', 'eval "$eec_runtime"',
                "cd {}".format(shlex.quote(str(snapshot))),
                "root-config --version", "make -j{}".format(jobs)]
    with log_path.open("w") as log:
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", "\n".join(commands)],
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("Payload build failed (exit {}); inspect {}".format(result.returncode, log_path))
    executable = snapshot / "eecAnalysis"
    if not executable.is_file() or not os.access(str(executable), os.X_OK):
        raise RuntimeError("Payload build did not produce an executable: {}".format(executable))
    for path in snapshot.rglob("*.o"):
        path.unlink()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=str(Path(__file__).absolute().parents[1]))
    parser.add_argument("--sample", choices=tuple(SAMPLES), default="OOMC")
    parser.add_argument("--input-dir")
    parser.add_argument("--submit-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--cmssw-base", default=CMSSW)
    parser.add_argument("--scram-arch", default=ARCH)
    parser.add_argument("--tag", default=timestamp().strftime("%Y%m%d_%H%M%S_%Z"))
    parser.add_argument("--background-methods", type=int, choices=(0, 3), default=0,
                        help="Regular modes: 0 = none, 3 = reflected and perpendicular cones; Response always uses 0")
    parser.add_argument("--match-eec-binning", action="store_true",
                        help="Use the EEC jet and track pT edges for sparse pT axes in every mode")
    parser.add_argument("--modes", nargs="+", choices=tuple(MODES), default=list(MODES))
    parser.add_argument("--build-jobs", type=int, default=2)
    parser.add_argument("--accounting-group", help="Optional explicitly chosen CERN accounting group")
    parser.add_argument("--list-only", action="store_true", help="Print input and job counts without creating files")
    args = parser.parse_args()
    args.input_dir = args.input_dir or "/eos/user/j/jvelazqu/Foresting/" + args.sample
    args.submit_dir = args.submit_dir or "/afs/cern.ch/user/j/jvelazqu/condor/production/{}/{}".format(args.sample, args.tag)
    args.output_dir = args.output_dir or "/eos/user/j/jvelazqu/HIN/eec-outputs/" + args.sample
    return args


def main():
    args = parse_args()
    repo = absolute_path(args.repo)
    input_dir = absolute_path(args.input_dir)
    submit_dir = absolute_path(args.submit_dir)
    output_dir = absolute_path(args.output_dir)
    cmssw_base = absolute_path(args.cmssw_base)
    modes = list(dict.fromkeys(args.modes))
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.tag):
        raise ValueError("--tag must contain only letters, numbers, underscores, dots, and hyphens")
    if not re.fullmatch(r"[A-Za-z0-9_]+", args.scram_arch):
        raise ValueError("Invalid SCRAM architecture")
    if args.build_jobs < 1:
        raise ValueError("--build-jobs must be positive")
    if args.accounting_group and not re.fullmatch(r"[A-Za-z0-9_.-]+", args.accounting_group):
        raise ValueError("Invalid accounting group")
    if not str(submit_dir).startswith("/afs/"):
        raise ValueError("--submit-dir must be under /afs; CERN submission files and logs must be on AFS")
    run_output = output_dir / args.tag
    eos_url(run_output)
    inputs = discover_inputs(input_dir, args.sample)
    job_count = len(inputs) * len(modes)
    summary = {"sample": args.sample, "inputs": len(inputs), "jobs": job_count, "modes": modes,
               "tag": args.tag, "submit_dir": str(submit_dir), "output_dir": str(run_output)}
    if args.list_only:
        print(json.dumps(summary, indent=2))
        return

    for name in ("jobs.sub", "manifest.json", "jobs.tsv", "analysis_payload.tar.gz"):
        if (submit_dir / name).exists():
            raise ValueError("Refusing to overwrite {}; choose a new --submit-dir".format(submit_dir / name))
    if run_output.exists():
        raise ValueError("Output campaign already exists: {}; choose a new --tag".format(run_output))
    if not cmssw_base.is_dir():
        raise ValueError("CMSSW release directory is unavailable: {}".format(cmssw_base))
    required = [repo / name for name in ("Makefile", "eecAnalysis.cxx", "src", "jetEnergyCorrections",
                                        "trackCorrectionTables", SAMPLES[args.sample]["card"],
                                        "condor/runOOMCJob.sh", "condor/checkOOMCOutput.C")]
    for path in required:
        if not path.exists():
            raise ValueError("Required repository file or directory is missing: {}".format(path))
    commit = git_output(repo, "rev-parse", "HEAD")
    git_status = git_output(repo, "status", "--porcelain")
    dirty = bool(git_status)
    provenance = commit + ("-dirty" if dirty else "")
    base_card = (repo / SAMPLES[args.sample]["card"]).read_text()
    cards = {mode: card_text(base_card, mode, args.background_methods, args.sample,
                             args.match_eec_binning) for mode in modes}
    submit_dir.mkdir(parents=True, exist_ok=True)
    (submit_dir / "logs").mkdir(exist_ok=True)
    (submit_dir / "cards").mkdir(exist_ok=True)
    print("Building the analysis once in local scratch; log: {}".format(submit_dir / "build.log"), flush=True)
    with tempfile.TemporaryDirectory(prefix="eec-{}-build-".format(args.sample.lower()), dir="/tmp") as scratch:
        snapshot = Path(scratch) / "analysis"
        snapshot.mkdir()
        for name in ("Makefile", "eecAnalysis.cxx"):
            shutil.copy2(str(repo / name), str(snapshot / name))
        for name in ("src", "jetEnergyCorrections", "trackCorrectionTables"):
            shutil.copytree(str(repo / name), str(snapshot / name), ignore=ignore_build_products)
        (snapshot / "condor").mkdir()
        shutil.copy2(str(repo / "condor/checkOOMCOutput.C"), str(snapshot / "condor/checkOOMCOutput.C"))
        (snapshot / "cards").mkdir()
        for mode, content in cards.items():
            (snapshot / "cards" / (mode + ".input")).write_text(content)
            (submit_dir / "cards" / (mode + ".input")).write_text(content)
        main_source = snapshot / "eecAnalysis.cxx"
        main_source.write_text(main_source.read_text().replace("GITHASHHERE", provenance))
        source_files = [{"path": path.relative_to(snapshot).as_posix(), "sha256": sha256(path),
                         "size_bytes": path.stat().st_size}
                        for path in sorted(snapshot.rglob("*")) if path.is_file()]
        source_digest = hashlib.sha256(json.dumps(source_files, sort_keys=True).encode()).hexdigest()
        build_snapshot(snapshot, cmssw_base, args.scram_arch, args.build_jobs, submit_dir / "build.log")
        executable_sha256 = sha256(snapshot / "eecAnalysis")
        archive = submit_dir / "analysis_payload.tar.gz"
        with tarfile.open(str(archive), "w:gz") as stream:
            stream.add(str(snapshot), arcname="analysis")
    archive_sha256 = sha256(archive)
    (submit_dir / "analysis_payload.tar.gz.sha256").write_text(archive_sha256 + "  analysis_payload.tar.gz\n")
    shutil.copy2(str(repo / "condor/runOOMCJob.sh"), str(submit_dir / "run_job.sh"))
    (submit_dir / "run_job.sh").chmod(0o755)
    rows = []
    for mode in modes:
        mode_output = run_output / mode
        mode_output.mkdir(parents=True, exist_ok=False)
        memory = 8000 if mode == "Response" else 4000
        for item in inputs:
            output_url = eos_url(mode_output / (item["stem"] + ".root"))
            rows.append("\t".join((mode, item["url"], output_url, item["stem"], str(memory))))
    (submit_dir / "jobs.tsv").write_text("\n".join(rows) + "\n")
    (submit_dir / "inputs.txt").write_text("\n".join(item["url"] for item in inputs) + "\n")
    accounting = '+AccountingGroup = "{}"\n'.format(args.accounting_group) if args.accounting_group else ""
    submit_text = """# Prepared {sample} campaign: {tag}; {count} inputs, {jobs} jobs.
# ROOT outputs are validated and copied by run_job.sh with XRootD.
universe = vanilla
initialdir = {workspace}
executable = {workspace}/run_job.sh
arguments = $(mode) $(input_url) $(output_url) $(stem) {cmssw} {arch} {checksum} {data_type}
log = {workspace}/logs/jobs.$(ClusterId).log
output = {workspace}/logs/$(mode).$(ClusterId).$(ProcId).out
error = {workspace}/logs/$(mode).$(ClusterId).$(ProcId).err
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {workspace}/analysis_payload.tar.gz
transfer_output_files = ""
request_cpus = 1
request_memory = $(memory_mb)
request_disk = 6000MB
requirements = (OpSysAndVer =?= "AlmaLinux9") && (Arch =?= "X86_64")
+JobFlavour = "workday"
+JobBatchName = "{sample}_{tag}"
{accounting}max_retries = 2
queue mode,input_url,output_url,stem,memory_mb from {workspace}/jobs.tsv
""".format(tag=args.tag, count=len(inputs), jobs=job_count, workspace=submit_dir, cmssw=cmssw_base,
           arch=args.scram_arch, checksum=archive_sha256, accounting=accounting,
           sample=args.sample, data_type=SAMPLES[args.sample]["data_type"])
    (submit_dir / "jobs.sub").write_text(submit_text)
    write_executable(submit_dir / "submit.sh", "#!/bin/bash\nset -euo pipefail\ncd {}\nexec condor_submit jobs.sub\n".format(shlex.quote(str(submit_dir))))
    command = ["python3", str(repo / "condor/prepareOOMC.py"), "--sample", args.sample, "--repo", str(repo),
               "--input-dir", str(input_dir), "--output-dir", str(output_dir),
               "--cmssw-base", str(cmssw_base), "--scram-arch", args.scram_arch,
               "--background-methods", str(args.background_methods), "--build-jobs", str(args.build_jobs),
               "--modes", *modes]
    if args.accounting_group:
        command.extend(["--accounting-group", args.accounting_group])
    if args.match_eec_binning:
        command.append("--match-eec-binning")
    prepare_script = """#!/bin/bash
set -euo pipefail
# Prepare production/{sample}/<tag>; a timestamped tag is generated unless supplied.
# Usage: ./prepare.sh [--tag NEW_TAG] [--submit-dir /afs/.../custom_workspace]
exec {command} "$@"
""".format(command=" ".join(shlex.quote(value) for value in command), sample=args.sample)
    write_executable(submit_dir / "prepare.sh", prepare_script)
    manifest = {"schema_version": 1, "created_at": timestamp().isoformat(), "prepared_only": True,
                "run_config": {"sample": args.sample, "data_type": SAMPLES[args.sample]["data_type"],
                               "base_card": SAMPLES[args.sample]["card"], "repo": str(repo), "input_dir": str(input_dir), "submit_dir": str(submit_dir),
                               "output_base": str(output_dir), "output_dir": str(run_output), "tag": args.tag,
                               "modes": modes, "background_methods": args.background_methods,
                               "match_eec_binning": args.match_eec_binning,
                               "cmssw_base": str(cmssw_base), "scram_arch": args.scram_arch,
                               "build_jobs": args.build_jobs, "accounting_group": args.accounting_group,
                               "job_flavour": "workday", "regular_memory_mb": 4000,
                               "response_memory_mb": 8000, "request_disk_mb": 6000},
                "input_count": len(inputs), "job_count": job_count,
                "inputs": inputs, "job_fields": ["mode", "input_url", "output_url", "stem", "memory_mb"],
                "provenance": {"git_commit": commit, "git_dirty": dirty,
                               "git_status_porcelain": git_status.splitlines(), "compiled_git_hash": provenance,
                               "source_manifest_sha256": source_digest, "source_files": source_files,
                               "executable_sha256": executable_sha256, "archive_sha256": archive_sha256,
                               "worker_sha256": sha256(submit_dir / "run_job.sh")}}
    (submit_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    readme = repo / "condor/README.md"
    if readme.is_file():
        shutil.copy2(str(readme), str(submit_dir / "README.md"))
    print(json.dumps(summary, indent=2))
    print("Prepared only. Review cards and manifest, then submit from AFS with:")
    print("  cd {} && ./submit.sh".format(shlex.quote(str(submit_dir))))
    print("Validate submission syntax without submitting:")
    print("  condor_submit -dry-run {} {}".format(shlex.quote(str(submit_dir / "dryrun.ads")),
                                                     shlex.quote(str(submit_dir / "jobs.sub"))))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(1)
