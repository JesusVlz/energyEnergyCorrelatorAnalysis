#!/bin/bash
set -Eeuo pipefail
trap 'eec_status=$?; printf "CAMPAIGN_FAILED line=%s status=%s\n" "$LINENO" "$eec_status" >&2; exit "$eec_status"' ERR
[[ $# -eq 4 ]]
eec_stage=$1
eec_index=$2
eec_system=$3
eec_sha=$4
[[ "$eec_sha" =~ ^[a-f0-9]{64}$ ]]
eec_scratch=$(cd "${_CONDOR_SCRATCH_DIR:-$PWD}" && pwd -P)
case "$eec_scratch" in /eos/*|/afs/*) echo 'Worker requires local scratch' >&2; exit 2;; esac
printf '%s  %s\n' "$eec_sha" "$eec_scratch/analysis_payload.tar.gz" | sha256sum -c -
eec_work=$(mktemp -d "$eec_scratch/run3campaign.XXXXXXXX")
tar -xzf "$eec_scratch/analysis_payload.tar.gz" -C "$eec_work"
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cd /cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_0_11
eec_runtime=$(scram runtime -sh)
eval "$eec_runtime"
set -u
cd "$eec_work/analysis"
ldd ./eecAnalysis > "$eec_work/libraries.txt"
if grep -Fq 'not found' "$eec_work/libraries.txt"; then cat "$eec_work/libraries.txt" >&2; exit 1; fi
printf 'RUN3_CAMPAIGN stage=%s system=%s index=%s ROOT=%s\n' "$eec_stage" "$eec_system" "$eec_index" "$(root-config --version)"
python3 validation/run3Campaign/batch_stage.py "$eec_stage" "$eec_index" "$eec_system" --scratch "$eec_work/work"
