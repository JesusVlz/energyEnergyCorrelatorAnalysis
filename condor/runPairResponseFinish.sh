#!/bin/bash
# Plain merge worker, following the same established payload and transfer pattern.
set -Eeuo pipefail
trap 'eec_status=$?; printf "PAIR_RESPONSE_FINISH_FAILED line=%s status=%s\n" "$LINENO" "$eec_status" >&2; exit "$eec_status"' ERR
[[ $# -eq 6 ]]
eec_kind=$1
eec_system=$2
eec_output_url=$3
eec_cmssw_base=$4
eec_arch=$5
eec_archive_sha=$6
case "$eec_kind" in pairs|response) ;; *) exit 2;; esac
case "$eec_system" in OO|pp) ;; *) exit 2;; esac
[[ "$eec_archive_sha" =~ ^[a-f0-9]{64}$ && "$eec_output_url" == root://* ]]
eec_scratch=$(cd "${_CONDOR_SCRATCH_DIR:-$PWD}" && pwd -P)
case "$eec_scratch" in /eos/*|/afs/*) exit 2;; esac
printf '%s  %s\n' "$eec_archive_sha" "$eec_scratch/analysis_payload.tar.gz" | sha256sum -c -
eec_work=$(mktemp -d "$eec_scratch/pair_response_finish.XXXXXXXX")
tar -xzf "$eec_scratch/analysis_payload.tar.gz" -C "$eec_work"
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH="$eec_arch"
cd "$eec_cmssw_base"
eec_runtime=$(scram runtime -sh)
eval "$eec_runtime"
set -u
cd "$eec_work/analysis"
python3 condor/finalize_pair_response.py --manifest "$PWD/campaign_manifest.json" --system "$eec_system" \
    --kind "$eec_kind" --work-dir "$eec_work/work" --output-dir "$eec_work/final"
tar -czf "$eec_work/${eec_kind}_final.tar.gz" -C "$eec_work/final" .
xrdcp --retry 3 --cksum adler32 --posc --force "$eec_work/${eec_kind}_final.tar.gz" "$eec_output_url"
printf 'PAIR_RESPONSE_FINISH_SUCCESS kind=%s system=%s output=%s\n' "$eec_kind" "$eec_system" "$eec_output_url"
