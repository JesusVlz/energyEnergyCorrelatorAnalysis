#!/bin/bash
set -uo pipefail

# --- 0. Entorno desde CVMFS (no depender del alias cmsenv ni de AFS) ---
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/
eval "$(scramv1 runtime -sh)"

# --- 1. Rutas ---
EOS_BASE="/eos/cms/store/group/phys_heavyions/jvelazqu/energycorr/eecAnalysis_OO_2025-09-10/eecOOHistograms/crab_eecAnalysis_OO_2025-09-10/260910_033733"
OUTPUT_DIR="/eos/cms/store/group/phys_heavyions/jvelazqu/energycorr/merged/"
OUTPUT_FILE="${OUTPUT_DIR}/merged_eecOO_2025-09-10.root"
mkdir -p "${OUTPUT_DIR}"

GOOD_LIST=$(mktemp)
BAD_LIST="${OUTPUT_DIR}/bad_files.txt"
: > "${BAD_LIST}"

# --- 2. Encontrar TODOS los archivos de output (find, no glob de shell) ---
echo "==> Buscando archivos en ${EOS_BASE}/000*"
mapfile -t ALL_FILES < <(find "${EOS_BASE}"/000* -name 'eecAnalysis_OO_2025-09-10_*.root' 2>/dev/null | sort -V)
echo "==> Encontrados ${#ALL_FILES[@]} archivos candidatos"

if [ "${#ALL_FILES[@]}" -eq 0 ]; then
    echo "ERROR: no se encontró ningún archivo. ¿Ruta o sufijo correctos?"
    exit 1
fi

# --- 3. Validar CADA archivo con ROOT antes de mergear ---
echo "==> Validando integridad (zombies y TTree truncados)..."
N_OK=0; N_BAD=0
for f in "${ALL_FILES[@]}"; do
    # 0 bytes: descarte inmediato sin abrir
    if [ ! -s "${f}" ]; then
        echo "${f}" >> "${BAD_LIST}"; N_BAD=$((N_BAD+1)); continue
    fi
    # zombie o recuperado (job muerto a media escritura)
    valid=$(root -l -b -q -e "
        TFile* f=TFile::Open(\"${f}\");
        if(!f||f->IsZombie()||f->TestBit(TFile::kRecovered)){printf(\"BAD\n\");}
        else{printf(\"OK\n\"); f->Close();}" 2>/dev/null | grep -E 'OK|BAD' | tail -1)
    if [ "${valid}" = "OK" ]; then
        echo "${f}" >> "${GOOD_LIST}"; N_OK=$((N_OK+1))
    else
        echo "${f}" >> "${BAD_LIST}"; N_BAD=$((N_BAD+1))
    fi
done

echo "==> Válidos: ${N_OK} | Malos: ${N_BAD}"
if [ "${N_BAD}" -gt 0 ]; then
    echo "==> ATENCIÓN: ${N_BAD} archivos inválidos listados en ${BAD_LIST}"
    echo "==> Esos eventos NO entran al merge. Reprocésalos antes de usar el resultado como física."
fi

# --- 4. Merge SOLO desde la lista validada (@file lee rutas de un fichero) ---
echo "==> Mergeando ${N_OK} archivos válidos..."
hadd -f -j 10 "${OUTPUT_FILE}" "@${GOOD_LIST}"
RC=$?
rm -f "${GOOD_LIST}"

# --- 5. Verificación ---
if [ ${RC} -eq 0 ]; then
    echo "========================================================="
    echo "Merge completado: ${OUTPUT_FILE}"
    [ "${N_BAD}" -gt 0 ] && echo "OJO: faltan ${N_BAD} archivos (ver ${BAD_LIST})"
    echo "========================================================="
else
    echo "ERROR: hadd falló con código ${RC}."
    exit 1
fi