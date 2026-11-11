#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 || $2 != card=* || $3 != output=* || $4 != location=* || $5 != mixinglist=* ]]; then
  echo "Usage: $0 JOB_ID card=FILE output=FILE location=INDEX mixinglist=INDEX" >&2
  exit 64
fi
JOBID=$1
CARD=${2#card=}
OUTPUT=${3#output=}
LOCATION=${4#location=}
MIXINGLIST=${5#mixinglist=}
if [[ ! $JOBID =~ ^[1-9][0-9]*$ || ! $LOCATION =~ ^[0-4]$ || -z $OUTPUT || $OUTPUT == */* ]]; then
  echo "Invalid CRAB job ID, file location, or output basename" >&2
  exit 64
fi
# Older configurations used 'none'; atoi("none") previously selected index 0.
if [[ $MIXINGLIST == none ]]; then MIXINGLIST=0; fi
if [[ ! $MIXINGLIST =~ ^-?[0-9]+$ ]]; then
  echo "mixinglist must be a numeric legacy index (or 'none')" >&2
  exit 64
fi

JOB_DIRECTORY=$PWD
if [[ $CARD != /* ]]; then CARD="$JOB_DIRECTORY/$CARD"; fi
if [[ ! -s $CARD ]]; then
  echo "Missing or empty submitted card: $CARD" >&2
  exit 66
fi
CARD_NAME=$(basename "$CARD")
WORKER_DIRECTORY=$(mktemp -d "$JOB_DIRECTORY/eec-worker.XXXXXX")
echo "CRAB job $JOBID: building in $WORKER_DIRECTORY"

# Isolate archive contents, including any historical cards or build artifacts.
tar -xzf "$JOB_DIRECTORY/eecOO.tar.gz" -C "$WORKER_DIRECTORY"
INPUT_DIRECTORY=$(mktemp -d "$WORKER_DIRECTORY/crab-inputs.XXXXXX")
tar -xf "$JOB_DIRECTORY/input_files.tar.gz" -C "$INPUT_DIRECTORY"
cp "$INPUT_DIRECTORY/job_input_file_list_${JOBID}.txt" "$WORKER_DIRECTORY/"
# The separately submitted card is authoritative, even if the archive has a copy.
cp "$CARD" "$WORKER_DIRECTORY/$CARD_NAME"

make -C "$WORKER_DIRECTORY" clean
make -C "$WORKER_DIRECTORY" eecAnalysis
cd "$WORKER_DIRECTORY"
./eecAnalysis "$JOBID" "$CARD_NAME" "$WORKER_DIRECTORY/$OUTPUT" "$LOCATION" false "$MIXINGLIST"
if [[ ! -s $OUTPUT ]]; then
  echo "Analyzer did not produce a nonempty output: $OUTPUT" >&2
  exit 65
fi
# Stage only a successful result into CRAB's original output-collection directory.
mv "$OUTPUT" "$JOB_DIRECTORY/$OUTPUT"
