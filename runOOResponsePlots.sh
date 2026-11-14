#!/usr/bin/env bash
# Sequential C++/ROOT response diagnostics from the already-projected campaign.
set -euo pipefail
PROJ="${PROJ:-${PWD}/projected/Run3MC_R4_20261112_v1_complete/OOMC}"
OUT="${OUT:-${PWD}/figures/OO_MC_responses_$(date +%Y%m%d_%H%M%S)}"
[[ $# -gt 0 ]] || set -- pairs jets unfoldingInputs
command -v root >/dev/null || { echo 'Initialize ROOT/CMSSW first.' >&2; exit 1; }
[[ -f plotting/Run3ResponsePlotTools.h ]] || { echo 'Run from the updated repository root.' >&2; exit 1; }
# Validate every input before creating output or running jobs.
for family in "$@"; do
  case "$family" in
    pairs) cfg=trackPairResolution; powers='1 2';;
    jets) cfg=jetResponse; powers='1';;
    unfoldingInputs) cfg=truthGuard5; powers='1 2';;
    *) echo "Unknown family: $family (pairs, jets, unfoldingInputs)" >&2; exit 1;;
  esac
  for scope in integrated split; do
    name="OO_MC_${scope}_${cfg}"
    for power in $powers; do
      suffix="_n${power}_processed.root"; [[ "$family" != jets ]] || suffix='_processed.root'
      file="${PROJ}/${name}/${name}${suffix}"
      [[ -r "$file" && -s "$file" ]] || { echo "Missing or empty input: $file" >&2; exit 1; }
    done
  done
done
[[ ! -e "$OUT" ]] || { echo "Choose a fresh OUT directory: $OUT" >&2; exit 1; }
mkdir -p "$OUT/logs"
make -C plotting 2>&1 | tee "$OUT/logs/build.log"
{ git rev-parse HEAD; root-config --version; git status --short; printf 'PROJ=%s\nOUT=%s\n' "$PROJ" "$OUT"; sha256sum plotting/*Run3.C plotting/Run3ResponsePlotTools.h; } > "$OUT/environment.txt"
for family in "$@"; do
  case "$family" in
    pairs) cfg=trackPairResolution; macro=trackParticleMatchingPlotterRun3; powers='1 2';;
    jets) cfg=jetResponse; macro=plotJetPtResponseMatrixRun3; powers='1';;
    unfoldingInputs) cfg=truthGuard5; macro=plotUnfoldingInputsRun3; powers='1 2';;
  esac
  for scope in integrated split; do
    name="OO_MC_${scope}_${cfg}"
    ranges=('-1 100'); [[ "$scope" != split ]] || ranges=('0 30' '30 50' '50 90')
    for power in $powers; do
      suffix="_n${power}_processed.root"; [[ "$family" != jets ]] || suffix='_processed.root'
      file="${PROJ}/${name}/${name}${suffix}"
      for range in "${ranges[@]}"; do
        read -r low high <<< "$range"
        target="$OUT/$family/$scope/C${low}-${high}/n${power}"
        log="$OUT/logs/${family}_${scope}_C${low}-${high}_n${power}.log"
        root -l -b -q "plotting/${macro}.C(\"${file}\",\"${target}\",${low}.,${high}.)" 2>&1 | tee "$log"
        grep -q RUN3_RESPONSE_PLOTTING_PASS "$log" || { echo "No success marker: $log" >&2; exit 1; }
      done
    done
  done
done
printf '\nFinished diagnostic plotting, not physics qualification.\nOutput: %s\n' "$OUT"
