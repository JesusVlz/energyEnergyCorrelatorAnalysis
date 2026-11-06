#!/bin/bash
# Same payload/runtime/input/output pattern as the working OOMC and ppMC workers.
set -Eeuo pipefail
trap 'eec_status=$?; printf "PAIR_RESPONSE_JOB_FAILED line=%s status=%s\n" "$LINENO" "$eec_status" >&2; exit "$eec_status"' ERR
[[ $# -eq 9 ]]
eec_mode=$1
eec_input_url=$2
eec_output_url=$3
eec_stem=$4
eec_cmssw_base=$5
eec_arch=$6
eec_archive_sha=$7
eec_data_type=$8
eec_task=$9
[[ "$eec_mode" == PairResponse ]]
case "$eec_data_type" in 10|12) ;; *) exit 2;; esac
[[ "$eec_stem" =~ ^(OO|pp)_[0-9]+$ ]]
[[ "$eec_archive_sha" =~ ^[a-f0-9]{64}$ && "$eec_task" =~ ^[0-9]+$ ]]
[[ "$eec_input_url" == root://* && "$eec_output_url" == root://* ]]
eec_scratch=$(cd "${_CONDOR_SCRATCH_DIR:-$PWD}" && pwd -P)
case "$eec_scratch" in /eos/*|/afs/*) echo 'Run workers in local scratch; submit from AFS.' >&2; exit 2;; esac
printf '%s  %s\n' "$eec_archive_sha" "$eec_scratch/analysis_payload.tar.gz" | sha256sum -c -
eec_work=$(mktemp -d "$eec_scratch/pair_response.XXXXXXXX")
tar -xzf "$eec_scratch/analysis_payload.tar.gz" -C "$eec_work"
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH="$eec_arch"
cd "$eec_cmssw_base"
eec_runtime=$(scram runtime -sh)
eval "$eec_runtime"
set -u
cd "$eec_work/analysis"
ldd ./eecAnalysis > "$eec_work/runtime-libraries.txt"
if grep -Fq 'not found' "$eec_work/runtime-libraries.txt"; then cat "$eec_work/runtime-libraries.txt" >&2; exit 1; fi
printf 'PAIR_RESPONSE_JOB input=%s stem=%s ROOT=%s\n' "$eec_input_url" "$eec_stem" "$(root-config --version)"
if [[ "${EEC_PREFLIGHT_ONLY:-0}" == 1 ]]; then echo 'PREFLIGHT_PASS'; exit 0; fi
python3 condor/pair_response_worker.py --analysis-dir "$PWD" --manifest "$PWD/campaign_manifest.json" \
    --task-index "$eec_task" --input-url "$eec_input_url" --input-local "$eec_work/input.root" --output-dir "$eec_work/products" --scratch "$eec_work/native"
# One completed artifact per forest, using the established direct POSC copy.
# Its two receipts bind every pair/response file; no remote directory creation.
[[ -f "$eec_work/products/pairs/stage_complete.json" && -f "$eec_work/products/response/stage_complete.json" ]]
tar -czf "$eec_work/$eec_stem.tar.gz" -C "$eec_work/products" pairs response
xrdcp --retry 3 --cksum adler32 --posc --force "$eec_work/$eec_stem.tar.gz" "$eec_output_url"
printf 'PAIR_RESPONSE_JOB_SUCCESS output=%s stem=%s\n' "$eec_output_url" "$eec_stem"
