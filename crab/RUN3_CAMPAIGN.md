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
an `allCones.input.in` card with an independently checked donor list. Do not
submit the `.input.in` template. The OO 35-file pilot pool is MB. The 35
HardProbes files in `ppRef2024_pool_pilot.txt` remain reserved out of the
versioned signal partitions; they are not the new ppRef mixing input. The
owner-selected five ZeroBiasPlusForward1 forests in
`ppRef2024_ZeroBiasPlusForward1_Nebraska_v2.txt` are a separate Nebraska
donor candidate. CRAB jobs still run at Vanderbilt and read those donors
remotely. Only mixed-cone cards check their selected donor list against the
signal files; event-level independence and tree contents must be checked
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

The ready OO all-cone card binds the existing 35-file MB pool reserved by
`splitOOForestFiles.py`, at
`cards/R4/materialized/OO_Data_split_allCones_MB35.input`. The original ten
round-robin parts and their two chunks remain unchanged: 20 tasks, at most
10,000 signal files per task, one file per job. The ppRef four-part, one-chunk
partition also remains unchanged. For both complete data campaigns, run the
plan commands first and then their submission counterparts:

```bash
python3 crab/submitRun3.py oo-data \
  --card cards/R4/materialized/OO_Data_split_allCones_MB35.input
python3 crab/submitRun3.py ppref-data \
  --card cards/R4/materialized/ppRef_Data_inclusive_allCones_NebraskaZB_v2.input
python3 crab/submitRun3.py oo-data \
  --card cards/R4/materialized/OO_Data_split_allCones_MB35.input --submit
python3 crab/submitRun3.py ppref-data \
  --card cards/R4/materialized/ppRef_Data_inclusive_allCones_NebraskaZB_v2.input --submit
```

The OO MB file list was selected for the pilot pool. Check its forest content
and per-job runtime with the pilot before interpreting full mixed-cone output.
The file-level overlap checks run for every part and chunk at submission.

The ready ppRef all-cone card binds the five owner-selected Nebraska ZeroBias
donors at `mixingFileList/ppRef2024_ZeroBiasPlusForward1_Nebraska_v2.txt`.
`BackgroundMethods 7` writes reflected, perpendicular, and mixed-cone
components in one analysis pass. It is a materialized card under
`cards/R4/materialized/ppRef_Data_inclusive_allCones_NebraskaZB_v2.input`.
The two full-data commands above use the existing Vanderbilt signal
partitions and the Nebraska ZeroBias donor URLs. Run them from lxplus with
`cmsenv`, a valid CMS proxy, and a rebuilt `crab/eecOO.tar.gz`. Each command
stops at the first failed CRAB submission; inspect the printed CRAB project
path before retrying. `--parts` can select remaining complete partitions.

The CRAB checks accept a nonempty donor list of any size and reject duplicate
or file-level overlapping donors. Check that all five Nebraska files open and
contain the needed trees and branches before adding `--submit`. The worker
archive must be rebuilt after adding or changing a donor list. If the local
`ppRef2024_pool_pilot.txt` was overwritten with the five new paths, compare
it with the separately tracked ZeroBias list before restoring the original
35-file reservation required by `splitPPRefForestFiles.py`.

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

## Complete configured MC set

The explicit driver below plans 20 OO MC and 13 ppRef MC tasks. It includes
OO integrated and split measurement, ppRef inclusive measurement, the
RecoGen/GenReco/GenGen closure cards, 5 GeV truth-guard response, jet and
track-pair QA, and Nominal/Loose/Tight RecoReco/RecoGen pair preparation.
The OO RecoReco measurement uses `BackgroundMethods 3` for both same-event
cones. The ppRef MC measurement uses `BackgroundMethods 0`: no mixed cone or
other background method is requested for ppRef MC. The existing
`mixingFileList/OO2025_HIJING_MC_pool_pilot_A.txt` references the pure-HIJING
donors, but its header requires an event-level partition and campaign-wide
reservation. The ordinary RecoReco `allCones` card cannot enforce that contract;
this CRAB driver rejects MC mixed-cone cards. The separate strict OO GenGen
mixed+perpendicular closure workflow (`BackgroundMethods 6`) uses the fixed
100,000-event HIJING donor prefix with `OOClosurePartitionPolicy 1` and writes
partner/exposure ledgers; see `condor/README.md` for its submitted campaigns.
The pair-preparation outputs do not turn on a residual pair correction.

```bash
python3 crab/buildEOSMCLists.py
python3 crab/submitRun3MCComplete.py
python3 crab/submitRun3MCComplete.py --submit
```

Run from the repository root with `cmsenv`, a valid CMS proxy, and a rebuilt
`crab/eecOO.tar.gz`. Submission stops at the first CRAB failure. To resume,
use `--sample` and `--stage` to select pending groups, or submit a pending
card individually with `crab/submitRun3.py --card`; existing CRAB projects
are not overwritten. The same full MC forests are processed once per selected
analysis card, so review the printed task count and CRAB quotas.

For each card, submit a small pilot with `--pilot 10` (or another 1..100
files). It writes a private list under `crab/inputLists/` with a distinct
`pilotN` task name. Add `--submit` only after reviewing the printed plan and the
card. The scripts do not check ROOT tree contents, remote file accessibility,
CRAB credentials, physics closure or per-job runtime; inspect a pilot output
and CRAB status before sending all data parts. The old Run 2 CRAB files live
under `crab/archive/run2/` for reference.
