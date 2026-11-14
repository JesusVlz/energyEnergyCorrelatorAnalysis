#!/usr/bin/env bash
# Run from energyEnergyCorrelatorAnalysis root after installing the new macro.
# This only READS the projected files. No reanalysis/reprojection is performed.
set -euo pipefail
PROJ="${PROJ:-${PWD}/projected/Run3MC_R4_20261112_v1_complete/OOMC}"
OUT="${OUT:-${PWD}/figures/OO_MC_tracking_$(date +%Y%m%d_%H%M%S)}"
if [ "$#" -eq 0 ]; then set -- Nominal; fi
[[ -f plotting/trackClosurePlotterRun3.C && -f plotting/Makefile ]] || {
  echo 'Run from the repository root and install plotting/trackClosurePlotterRun3.C first.' >&2; exit 1;
}
command -v root >/dev/null || { echo 'Initialize your ROOT/CMSSW environment first.' >&2; exit 1; }
for wp in "$@"; do
  case "$wp" in Nominal|Loose|Tight) ;; *) echo "Invalid working point: $wp" >&2; exit 1;; esac
  for mode in RecoReco RecoGen; do
    for power in 1 2; do
      cfg="OO_${wp}_${mode}"
      input="${PROJ}/${cfg}/${cfg}_n${power}_processed.root"
      [[ -r "$input" && -s "$input" ]] || { echo "Missing/empty input: $input" >&2; exit 1; }
    done
  done
done
# A fresh directory prevents confusion with previous plots and partial runs.
if [[ -e "$OUT" ]]; then echo "Output already exists; set OUT to a new directory: $OUT" >&2; exit 1; fi
mkdir -p "$OUT/logs"
make -C plotting 2>&1 | tee "$OUT/logs/build.log"
{
  printf 'Repository: '; pwd
  printf 'Commit: '; git rev-parse HEAD
  printf 'ROOT: '; root-config --version
  printf '\nLocal changes:\n'; git status --short
  printf '\nMacro hash:\n'; sha256sum plotting/trackClosurePlotterRun3.C
  printf '\nPROJ=%s\nOUT=%s\n' "$PROJ" "$OUT"
} > "$OUT/environment.txt"
for wp in "$@"; do
  rr="OO_${wp}_RecoReco"
  rg="OO_${wp}_RecoGen"
  for power in 1 2; do
    reco="${PROJ}/${rr}/${rr}_n${power}_processed.root"
    gen="${PROJ}/${rg}/${rg}_n${power}_processed.root"
    tracks=false; [[ "$power" -eq 1 ]] && tracks=true
    for range in '0 30' '30 50' '50 90'; do
      read -r low high <<< "$range"
      out="${OUT}/${wp}/C${low}-${high}/n${power}"
      log="${OUT}/logs/${wp}_C${low}-${high}_n${power}.log"
      root -l -b -q \
        "plotting/trackClosurePlotterRun3.C(\"${reco}\",\"${gen}\",\"${out}\",${low}.,${high}.,2.,\"${wp}\",${tracks},true)" \
        2>&1 | tee "$log"
      grep -q 'RUN3_TRACKING_PLOTTING_PASS' "$log" || {
        echo "No successful completion marker; inspect $log" >&2; exit 1;
      }
    done
  done
done
printf '\nFinished requested plotting jobs. Physics closure is not inferred from execution.\nOutput: %s\n' "$OUT"
