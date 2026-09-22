#!/bin/bash
set -uo pipefail

FILELIST="$1"   # UNA filelist, la que Condor asignó a este job

cd /afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/
cmsenv

OUTPUT_DIR="/eos/cms/store/group/phys_heavyions/jvelazqu/energycorr/condor_R8_axisWTA_inclusive"
mkdir -p "${OUTPUT_DIR}"

echo "==> Job procesando filelist: ${FILELIST}"

while IFS= read -r INPUT_FILE || [ -n "${INPUT_FILE}" ]; do
    [ -z "${INPUT_FILE}" ] && continue
    case "${INPUT_FILE}" in \#*) continue ;; esac

    FILE_NAME=$(basename "${INPUT_FILE}")
    BASE="${FILE_NAME%.*}"
    EXT="${FILE_NAME##*.}"
    OUTPUT_FILE="${OUTPUT_DIR}/${BASE}_processed.${EXT}"

    echo "----> ${FILE_NAME} -> ${BASE}_processed.${EXT}"
    ./eecAnalysis "${INPUT_FILE}" cardEECOO_inclusive.input "${OUTPUT_FILE}" 1 true
done < "${FILELIST}"

echo "==> Filelist completada: ${FILELIST}"
