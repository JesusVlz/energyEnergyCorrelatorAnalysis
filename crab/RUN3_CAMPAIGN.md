# Run 3 EEC CRAB campaign (R=0.4)

Run these commands in the `energyEnergyCorrelatorAnalysis` checkout on lxplus with
`cmsenv` and a valid CMS proxy. They prepare jobs; nothing submits without
`--submit`. The current R4 catalog is in `cards/R4/README.md`.

## Inputs and execution

| Sample | Input list | Default card | Tasks |
|---|---|---|---:|
| OO data | Vanderbilt `OOForestfiles/OOForestFiles_partNN_chunkNN.txt` | `cards/R4/measurement/OO_Data_split_sameEventCones.input` | 20 |
| ppRef data | Vanderbilt `ppRefForestfiles/ppRefForestFiles_partNN_chunk01.txt` | `cards/R4/measurement/ppRef_Data_inclusive_sameEventCones.input` | 4 |
| OO MC | generated `inputLists/OOMC_EOS.txt` | `cards/R4/measurement/OO_MC_split_none.input` | 1 |
| ppRef MC | generated `inputLists/ppRefMC_EOS.txt` | `cards/R4/measurement/ppRef_MC_inclusive_none.input` | 1 |

All tasks run **only at `T2_US_Vanderbilt`**. The data forests remain at
Vanderbilt; the MC inputs are read remotely through `root://eosuser.cern.ch`
from the user's EOS directories. CRAB stages histogram output at
`T2_CH_CERN` (`/store/group/phys_heavyions/jvelazqu/eec/`). The whitelist is
an execution-site restriction, not a stageout-site setting.
Versioned data lists may contain LFNs or XRootD PFNs; the data configuration
normalizes each to its `/store/` LFN, and the analyzer reads it from the
Vanderbilt XRootD endpoint (`location=2`). The MC PFNs are preserved.

Build input lists from the *actual* EOS directory contents rather than copying
the older MC lists that point at Vanderbilt. The OO glob is
`unembedded_L1HLT_objects_MC_*.root` under `/eos/user/j/jvelazqu/Foresting/OOMC`;
the ppRef glob is `HiForestMiniAOD_MC_*.root` under
`/eos/user/j/jvelazqu/Foresting/ppMC`. Check that all intended MC forest files
follow these patterns. The generated lists are kept local and ignored by Git.

```bash
cd /eos/home-j/jvelazqu/HIN/CMSSW_15_0_11/src/energyEnergyCorrelatorAnalysis
cmsenv
voms-proxy-init -voms cms
python3 crab/buildEOSMCLists.py
python3 cards/build_r4_catalog.py --check
python3 crab/splitOOForestFiles.py
python3 crab/splitPPRefForestFiles.py
python3 crab/makeOOAnalysisTar.py
```

If `--check` or a splitter fails, examine the checked-out revision and input
lists before submission. Use a checked-out commit and matching worker archive;
rebuild the archive after any code or calibration change.

## Plan, pilot, submit

The default data cards calculate the same-event reflected and perpendicular
cones together in one pass, without a donor pool. For mixed cones, materialize
an `allCones.input.in` card with an independently validated donor list. Do not
submit the `.input.in` template. The current OO 35-file pool is a pilot MB
sample; the ppRef pilot donor pool consists of HardProbes-filtered forests and
does not establish an unbiased ZeroBias reference. The config validates
disjoint signal and pool files, but event-level independence must be checked
separately. Retaining one data file per CRAB job isolates empty input forests;
each task stays at or below the CRAB 10,000-job limit. Do not run every
background mask across the full data sample.

When a suitable pool has been validated, make a runnable all-cone card under
`cards/R4/materialized/` and select it with `--card`. For the existing OO
pilot pool, the materialization commands are:

```bash
mkdir -p cards/R4/materialized
printf '{"MIXING_LIST":"mixingFileList/OO2025_MB_pool_pilot.txt"}\n' > cards/R4/materialized/oo-pool.json
python3 cards/build_r4_catalog.py materialize \
  cards/R4/measurement/OO_Data_split_allCones.input.in \
  --values cards/R4/materialized/oo-pool.json \
  --output cards/R4/materialized/OO_Data_split_allCones.input
python3 crab/submitRun3.py oo-data --parts 1 \
  --card cards/R4/materialized/OO_Data_split_allCones.input
```

The analogous ppRef pool must first pass the donor provenance checks stated
above. Do not substitute the HardProbes pilot list as a certified unbiased
ZeroBias donor sample. The card and its actual donor list must be included
when rebuilding the worker tarball.

```bash
python3 crab/submitRun3.py oo-data --pilot 10
python3 crab/submitRun3.py ppref-data --pilot 10
python3 crab/submitRun3.py oo-mc --pilot 10
python3 crab/submitRun3.py ppref-mc --pilot 10

# After a successful pilot, append --submit for the actual tasks:
python3 crab/submitRun3.py oo-data --pilot 10 --submit
python3 crab/submitRun3.py ppref-data --pilot 10 --submit
python3 crab/submitRun3.py oo-mc --pilot 10 --submit
python3 crab/submitRun3.py ppref-mc --pilot 10 --submit

# Once the pilot is checked, plan all production parts, then submit:
python3 crab/submitRun3.py oo-data
python3 crab/submitRun3.py ppref-data
python3 crab/submitRun3.py oo-data --submit
python3 crab/submitRun3.py ppref-data --submit
python3 crab/submitRun3.py oo-mc --submit
python3 crab/submitRun3.py ppref-mc --submit
```

An existing project is never silently reused. Set `--tag R4v2` or another
unique short tag for a replacement, and check its output before proceeding.
`--card` selects another runnable `.input` file under `cards/R4/`. The request
name contains its basename, so different card outputs do not overwrite one
another. For example, the complementary MC studies can be planned and sent
individually, using the same EOS file lists:

```bash
python3 crab/submitRun3.py oo-mc --card cards/R4/jet-qa/OO_MC_split_jetResponse.input
python3 crab/submitRun3.py oo-mc --card cards/R4/mc-closure/OO_MC_split_RecoGen.input
python3 crab/submitRun3.py oo-mc --card cards/R4/mc-closure/OO_MC_split_GenReco.input
python3 crab/submitRun3.py oo-mc --card cards/R4/mc-closure/OO_MC_split_GenGen.input
python3 crab/submitRun3.py oo-mc --card cards/R4/response/OO_MC_split_truthGuard5.input
python3 crab/submitRun3.py ppref-mc --card cards/R4/jet-qa/ppRef_MC_inclusive_jetResponse.input
python3 crab/submitRun3.py ppref-mc --card cards/R4/response/ppRef_MC_inclusive_truthGuard5.input
```

The nominal/loose/tight RecoReco and RecoGen pair preparation cards in
`cards/R4/pairs/` are also selectable. Run them as separate bounded studies;
the derived maps require training and holdout checks before a
`corrected-templates` card can be materialized. The response cards produce
ingredients for unfolding, not unfolded spectra. The `truthGuard10` cards
are sensitivity checks. Additional track-pair and matched-jet diagnostics
are available in the same R4 catalog. R8 and input-dependent systematic
templates remain outside this automatic campaign pending their inputs and
calibration.

For each card, submit a small pilot with `--pilot 10` (or another 1..100
files). It writes a private list under `crab/inputLists/` with a distinct
`pilotN` task name. Add `--submit` only after reviewing the printed plan and the
card. The scripts do not check ROOT tree contents, remote file accessibility,
CRAB credentials, physics closure or per-job runtime; inspect a pilot output
and CRAB status before sending all data parts. The old Run 2 CRAB files live
under `crab/archive/run2/` for reference.
