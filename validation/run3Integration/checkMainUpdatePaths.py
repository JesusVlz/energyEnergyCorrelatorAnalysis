#!/usr/bin/env python3
"""Prevent the one-time integration from committing unrelated changes."""
from pathlib import Path
import subprocess
root=Path(__file__).resolve().parents[2]
allowed={
    'src/EECAnalyzer.cxx','src/HighForestReader.cxx','src/HighForestReader.h',
    'src/UnfoldingForestReader.cxx','cards/build_r4_catalog.py','cards/README.md',
    'README.md','validation/run3Pairs/export_forest.py',
    'validation/run3Response/validateUnfoldingForestReader.cxx',
    'validation/run3Workflow/check_cards.py',
    '.github/workflows/oo-integration-validation.yml',
}
prefixes=('plotting/deriveMonteCarloWeights', 'plotting/testDeriveMonteCarloWeights',
          'plotting/README_Run3MCWeights.md','validation/run3Integration/','cards/R4/')
paths=subprocess.check_output(['git','diff','--name-only','HEAD'],cwd=root,text=True).splitlines()
paths+=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=root,text=True).splitlines()
for p in paths:
    if p in allowed or p.startswith(prefixes):continue
    if p.endswith('.input') and not any('/'+s+'/' in '/'+p for s in ('archive','evidence','results','pending')):continue
    raise RuntimeError('Unreviewed changed path: '+p)
print('PASS: reviewed integration paths only')
