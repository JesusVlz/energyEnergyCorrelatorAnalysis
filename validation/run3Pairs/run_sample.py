#!/usr/bin/env python3
"""Run a fixed local R4 MC sample through the native analyzer-backed exporter."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import subprocess
import time

REPO = Path(__file__).resolve().parents[2]
SAMPLES = {
    "OO": ("/eos/user/j/jvelazqu/Foresting/OOMC", "unembedded_L1HLT_objects_MC_", "OO2025_PHOO_pythia_embedded_2GeV_constit_cut_260926"),
    "pp": ("/eos/user/j/jvelazqu/Foresting/ppMC", "HiForestMiniAOD_MC_", "ppRef2024_Pythia_QCD_2GeV_constit_cut_260926"),
}


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8*1024*1024), b""):
            result.update(block)
    return result.hexdigest()



def reusable_result(task, source_sha256):
    """Reuse only a successful export bound to the same request and file bytes."""
    target = Path(task["output"])
    receipt = target.with_suffix(".run.json")
    metadata_path = Path(str(target)+".metadata.json")
    try:
        saved = json.loads(receipt.read_text())
        if (saved.get("exit_code") != 0 or saved.get("source_sha256") != source_sha256
                or any(saved.get(key) != value for key, value in task.items())
                or saved.get("input_sha256") != digest(Path(task["input"]))):
            return None
        metadata = json.loads(metadata_path.read_text())
        if metadata != saved.get("metadata") or metadata.get("output_sha256") != digest(target):
            return None
    except (OSError, ValueError):
        return None
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--oo-files", type=int, default=34)
    parser.add_argument("--pp-files", type=int, default=10)
    parser.add_argument("--events-per-file", type=int, default=30000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if min(args.oo_files,args.pp_files)<0 or not 1<=args.workers<=4 or args.events_per_file<1:
        parser.error("nonnegative file counts, positive event bound, 1..4 workers required")
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    tasks=[]
    for system, count in (("OO",args.oo_files),("pp",args.pp_files)):
        directory, prefix, namespace = SAMPLES[system]
        for number in range(1,count+1):
            source=Path(directory)/(prefix+str(number)+".root")
            if not source.is_file():
                raise FileNotFoundError(source)
            target=output/(f"{system}_{number:03d}.jsonl")
            command=["python3",str(REPO/"validation/run3Pairs/export_forest.py"),str(source),str(target),
                     "--system",system,"--sample-namespace",namespace,"--max-events",str(args.events_per_file)]
            tasks.append(dict(system=system, input=str(source),output=str(target),command=command))
    manifest={"procedure":"run3-an751-v1", "tasks":tasks,"events_per_file":args.events_per_file,
              "workers":args.workers,"training_validation":"unchanged namespaced physical-event FNV1a64 modulo2",
              "source_sha256":{str(p.relative_to(REPO)):digest(p) for p in [REPO/"eecAnalysis",REPO/"src/EECAnalyzer.cxx",REPO/"src/HighForestReader.cxx",REPO/"src/TrackPairEfficiencyCorrector.cxx",REPO/"validation/run3Pairs/export_forest.py"]}}
    (output/"sample_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    (output/"inputs.txt").write_text("\n".join(t["output"] for t in tasks)+"\n")
    for system in SAMPLES:
        (output/f"{system}_inputs.txt").write_text("\n".join(t["output"] for t in tasks if t["system"] == system)+"\n")
    def run(task):
        started=time.monotonic(); target=Path(task["output"]); receipt=target.with_suffix(".run.json")
        if args.resume:
            saved=reusable_result(task,manifest["source_sha256"])
            if saved is not None:
                return saved
        with target.with_suffix(".export.log").open("w") as log:
            code=subprocess.run(task["command"],cwd=REPO,stdout=log,stderr=subprocess.STDOUT).returncode
        result=dict(task,exit_code=code,elapsed_seconds=time.monotonic()-started,
                    input_sha256=digest(Path(task["input"])),source_sha256=manifest["source_sha256"])
        if code==0:
            result["metadata"]=json.loads(Path(str(target)+".metadata.json").read_text())
        receipt.write_text(json.dumps(result,indent=2)+"\n")
        print(json.dumps({k:result[k] for k in ("system","input","exit_code","elapsed_seconds")}),flush=True)
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results=list(pool.map(run,tasks))
    summary={"files":len(results),"failed":[r["input"] for r in results if r["exit_code"]],
             "source_sha256":manifest["source_sha256"],"systems":{}}
    for system in SAMPLES:
        counts={}
        for result in results:
            if result["system"]==system and result["exit_code"]==0:
                for key,value in result["metadata"]["counts"].items(): counts[key]=counts.get(key,0)+value
        summary["systems"][system]=counts
    (output/"sample_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary),flush=True)
    if summary["failed"]: raise SystemExit(1)


if __name__=="__main__": main()
