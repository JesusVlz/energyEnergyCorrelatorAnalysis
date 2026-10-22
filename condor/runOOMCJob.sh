#!/bin/bash
# Process one input forest and one O-O MC configuration in worker scratch.
set -Eeuo pipefail
trap 'eec_status=$?; printf "OOMC_JOB_FAILED line=%s status=%s\n" "$LINENO" "$eec_status" >&2; exit "$eec_status"' ERR

if [[ $# -lt 7 || $# -gt 8 ]]; then
    echo "Usage: $0 MODE INPUT_URL OUTPUT_URL STEM CMSSW_BASE SCRAM_ARCH ARCHIVE_SHA256 [EXPECTED_DATA_TYPE]" >&2
    exit 2
fi
eec_mode=$1
eec_input_url=$2
eec_output_url=$3
eec_stem=$4
eec_cmssw_base=$5
eec_arch=$6
eec_archive_sha=$7
eec_expected_data_type=${8:-10}

case "$eec_expected_data_type" in 10|12) ;; *)
    echo "Expected MC data type must be 10 (OO) or 12 (ppRef)." >&2
    exit 2
esac

case "$eec_mode" in RecoReco|RecoGen|GenReco|GenGen|Response) ;; *)
    echo "Unknown MC configuration: $eec_mode" >&2
    exit 2
esac
[[ "$eec_stem" =~ ^[A-Za-z0-9_.-]+$ && "$eec_stem" != "." && "$eec_stem" != ".." ]]
[[ "$eec_archive_sha" =~ ^[a-fA-F0-9]{64}$ ]]
[[ "$eec_arch" =~ ^[A-Za-z0-9_]+$ ]]
[[ "$eec_input_url" == root://* && "$eec_output_url" == root://* ]]

eec_scratch=$(cd "${_CONDOR_SCRATCH_DIR:-$PWD}" && pwd -P)
case "$eec_scratch" in /eos/*|/afs/*)
    echo "Run workers in local scratch; submit.sh belongs on AFS." >&2
    exit 2
esac
eec_archive="$eec_scratch/analysis_payload.tar.gz"
[[ -f "$eec_archive" ]]
printf '%s  %s\n' "$eec_archive_sha" "$eec_archive" | sha256sum -c -
eec_work=$(mktemp -d "$eec_scratch/oomc.XXXXXXXX")
tar -xzf "$eec_archive" -C "$eec_work"
[[ -x "$eec_work/analysis/eecAnalysis" ]]
[[ -f "$eec_work/analysis/cards/$eec_mode.input" ]]
[[ -f "$eec_work/analysis/condor/checkOOMCOutput.C" ]]

# cmsset_default and SCRAM runtime snippets can reference unset shell variables.
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH="$eec_arch"
cd "$eec_cmssw_base"
eec_runtime=$(scram runtime -sh)
eval "$eec_runtime"
set -u
cd "$eec_work/analysis"

ldd ./eecAnalysis > "$eec_work/runtime-libraries.txt"
if grep -Fq 'not found' "$eec_work/runtime-libraries.txt"; then
    cat "$eec_work/runtime-libraries.txt" >&2
    exit 1
fi
printf 'OOMC_JOB mode=%s host=%s ROOT=%s\n' "$eec_mode" "$(hostname)" "$(root-config --version)"
printf 'Input: %s\nOutput: %s\n' "$eec_input_url" "$eec_output_url"

if [[ "${EEC_PREFLIGHT_ONLY:-0}" == 1 ]]; then
    printf 'PREFLIGHT_PASS mode=%s payload=%s\n' "$eec_mode" "$eec_archive_sha"
    exit 0
fi

eec_input="$eec_work/input.root"
eec_output="$eec_work/$eec_stem.root"
xrdcp --retry 3 --cksum adler32 "$eec_input_url" "$eec_input"
./eecAnalysis "$eec_input" "cards/$eec_mode.input" "$eec_output" 1 true
root -l -b -q "condor/checkOOMCOutput.C(\"$eec_output\",\"$eec_mode\",$eec_expected_data_type)"

# Unique run/mode/input paths permit retries of the same job. POSC removes
# an incomplete destination if the transfer fails.
xrdcp --retry 3 --cksum adler32 --posc --force "$eec_output" "$eec_output_url"
printf 'OOMC_JOB_SUCCESS mode=%s output=%s\n' "$eec_mode" "$eec_output_url"
