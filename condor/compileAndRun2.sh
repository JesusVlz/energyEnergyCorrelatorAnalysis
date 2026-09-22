#!/bin/bash

INPUT_FILE=$1

# 2. Extraer el nombre base del archivo input
FILE_NAME=$(basename "${INPUT_FILE}")

# 2b. Separar nombre y extensión para insertar "processed"
BASE="${FILE_NAME%.*}"      # nombre sin extensión (ej: L1_object_data_575)
EXT="${FILE_NAME##*.}"      # extensión (ej: root)
OUT_NAME="${BASE}_processed.${EXT}"   # -> L1_object_data_575_processed.root

# 3. Definir el directorio de salida en EOS
OUTPUT_DIR="/eos/cms/store/group/phys_heavyions/jvelazqu/energycorr/condor_R4_axisWTA_cen/0006"
OUTPUT_FILE="${OUTPUT_DIR}/${OUT_NAME}"

# 4. Moverse al entorno CMSSW y cargar las librerías
cd /afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/

cmsenv

# 5. Ejecutar el análisis
./eecAnalysis "${INPUT_FILE}" cardEECOO.input "${OUTPUT_FILE}" 1 true

#./eecAnalysis /eos/cms/store/group/phys_heavyions/jvelazqu/OO2025_dataset_MinimumBias_filtered/IonPhysics0/OO_MinimumBias_filtered_IonPhysics0/260908_000643/0000/L1_object_data_541.root cardEECOO.input test_OO.root 1 true
