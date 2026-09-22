#!/bin/bash

# 1. Definir rutas de EOS
EOS_DIR="/eos/cms/store/group/phys_heavyions/jvelazqu/energycorr/condor_R4_axisWTA_inclusive"
OUTPUT_DIR="/eos/cms/store/group/phys_heavyions/jvelazqu/energycorr/merged/R4_axisWTA_inclusive"
OUTPUT_FILE="${OUTPUT_DIR}/merged_R4_axisWTA_inclusive.root"

# Crear la carpeta de salida si no existe
mkdir -p $OUTPUT_DIR

echo "==> Iniciando el proceso de unificación (Merge)..."
echo "==> Buscando archivos en: ${EOS_DIR}"

# 2. Cargar el entorno CMSSW/ROOT por si acaso no está activo en la terminal actual
cd /afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/
cmsenv

hadd -f -j 10 ${OUTPUT_FILE} ${EOS_DIR}/L1_object_data_*

# 5. Verificación final
if [ $? -eq 0 ]; then
    echo "========================================================="
    echo "¡Merge completado con éxito!"
    echo "Archivo final creado en: ${OUTPUT_FILE}"
    echo "========================================================="
else
    echo "Error: Ocurrió un problema durante la ejecución de hadd."
    exit 1
fi
