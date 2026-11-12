#!/usr/bin/env bash
# Software checks only. Does not submit production or claim physics closure.
set -euo pipefail
cd "$(dirname "$0")/../.."
EVIDENCE=${EEC_MAIN_EVIDENCE:-/tmp/run3-evidence}
mkdir -p "$EVIDENCE"
exec > >(tee "$EVIDENCE/checks.log") 2>&1
python validation/run3Integration/checkSlimTrackCards.py --write
python cards/build_r4_catalog.py --write
python cards/build_r4_catalog.py --check
python validation/run3Integration/checkSlimTrackCards.py
python validation/run3Workflow/check_cards.py
make -j2
make -C plotting -j2
CFLAGS=$(root-config --cflags)
LIBS=$(root-config --libs)
g++ $CFLAGS -fno-access-control -Isrc validation/run3Integration/testSlimTrackPolicy.cxx src/ForestReader.o src/HighForestReader.o src/UnfoldingForestReader.o $LIBS -o "$EVIDENCE/slim_tracks"
"$EVIDENCE/slim_tracks" "$EVIDENCE/slim-fixtures"
g++ $CFLAGS validation/run3Response/validateUnfoldingForestReader.cxx src/UnfoldingForestReader.o $LIBS -o "$EVIDENCE/unfolding_reader"
"$EVIDENCE/unfolding_reader" "$EVIDENCE/response-fixtures"
python -m unittest discover -s validation/run3Pairs -p 'test_*.py' -v
python validation/run3Pairs/prepare_inputs.py audit-tables --output "$EVIDENCE/tracking_tables.json"
python validation/run3Pairs/prepare_inputs.py cards --output "$EVIDENCE/pair-cards"
diff -r validation/run3Pairs/cards "$EVIDENCE/pair-cards"
python validation/run3Jets/run_selection_tests.py --output "$EVIDENCE/jets"
g++ -std=c++17 $CFLAGS validation/test_full_pf_multiplicities.cxx src/Run3JetSelection.cxx $LIBS -o "$EVIDENCE/pf_counts"
"$EVIDENCE/pf_counts"
python -m unittest discover -s validation/run3Unfolding -p 'test_study*.py' -v
mkdir -p "$EVIDENCE/aclic"
for macro in plotting/compareDifferentEECpairingsRun3.C plotting/validateBackgroundSubtractionRun3.C validation/run3Plotting/checkRun3Closure.C; do
  MACRO="$macro" BUILD="$EVIDENCE/aclic" root -l -b -q -e 'if(gSystem->Load("plotting/DrawingClasses.so")<0) gSystem->Exit(1); if(!gSystem->CompileMacro(gSystem->Getenv("MACRO"),"kf","",gSystem->Getenv("BUILD"))) gSystem->Exit(1);'
done
(cd plotting; root -l -b -q 'testDeriveMonteCarloWeightsRun3.C()')
cp -r /tmp/eec_weight_test_* "$EVIDENCE/" 2>/dev/null || true
python - <<'PY'
import ROOT
print('RooUnfold runtime available:',ROOT.gSystem.Load('libRooUnfold')>=0)
print('Software checks only; native MC and full-chain physics closure NOT_RUN.')
PY
git diff --check
