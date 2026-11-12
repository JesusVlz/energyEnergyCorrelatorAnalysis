# Running the code in crab

All the necessary configuration for running the code using the CRAB system is found from this folder. Below there is a simple example on how to run the analysis on CRAB and also detailed descriptions of all the files in this folder.

Active task lists, reserved donor pools and their manifests are versioned.
Worker tarballs, CRAB task directories, `submissionRecords/`, bytecode and logs
are local artifacts. Historical `previous_10parts/` directories mentioned below
are retained in the production workspace and excluded from Git. Rebuild the
worker from source with `python3 makeOOAnalysisTar.py` before submitting.

## OO mixed/perpendicular pilot: one signal file per job

`crabEECOO.py` uses `FileBased` splitting with `config.Data.unitsPerJob = 1`.
An empty signal forest can still fail its own job, but no healthy signal files
share that job. The analyzer and worker failure handling are unchanged.

The 132,863 signal files remain in their original ten disjoint 10% partitions,
with all 35 files in `../mixingFileList/OO2025_MB_pool_pilot.txt` excluded.
Each partition is split sequentially into two task lists: `chunk01` has 10,000
files and `chunk02` has 3,287 files for parts 01–03 or 3,286 for parts 04–10.
This gives 20 task lists in total. Both chunks of part01 together preserve the
original 13,287-file pilot sample.

The task size follows the 10,000-job default in the
[CRAB task splitter](https://github.com/dmwm/CRABServer/blob/8069e0cb72445e9a3a13dfb98469a0920210ed12/src/python/TaskWorker/Actions/Splitter.py#L23).
That limit is server-configurable; the splitter rejects tasks above the limit.
With one file per job, each generated list stays at or below 10,000 jobs.

Active lists are named `OOForestfiles/OOForestFiles_part01_chunk01.txt` through
`OOForestfiles/OOForestFiles_part10_chunk02.txt`. The original ten lists are
preserved under `OOForestfiles/previous_10parts/`. From `crab/`, run:

```bash
python3 splitOOForestFiles.py
```

The splitter writes the task lists and `OOForestfiles/manifest.json`.
An identical rerun verifies the existing outputs without rewriting them;
changed existing outputs are rejected. Configuration loading checks signal
duplicates and overlap with the reserved pool after converting XRootD URLs
to `/store/` LFNs.

The configuration defaults to part01, chunk01. Select a part (1–10) and chunk
(1–2) with `OO_CRAB_PART` and `OO_CRAB_CHUNK`. For example, these commands
submit the two chunks comprising the original 10% pilot:

```bash
OO_CRAB_PART=1 OO_CRAB_CHUNK=1 crab submit -c crabEECOO.py
OO_CRAB_PART=1 OO_CRAB_CHUNK=2 crab submit -c crabEECOO.py
```

Request and output names include the selected part and chunk, for example
`eecAnalysis_OO_mixedPerp_part01_chunk01_Run3_v2`, keeping replacement outputs
separate from the killed v1 task. No replacement tasks were submitted as part
of preparing this resplit.

The default execution site is `T2_US_Vanderbilt`. Set `additionalSites` in the
configuration to `['T2_US_Purdue', 'T2_US_Nebraska']` to permit those sites too.
A whitelist is an allowed-site list, not an ordered fallback preference.
Inputs still use Vanderbilt XRootD; outputs go to `T2_CH_CERN` under
`/store/group/phys_heavyions/jvelazqu/eec/`.
See the [CRAB configuration reference](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile).

The OO card selects mixed plus perpendicular cones (`BackgroundMethods 6`),
the explicit 35-file `MixingFileList`, regular forests (`MegaSkimMixing 0`),
and unbinned Vz/hiBin matching (`EventMatchForMixing 0`, `BinnedMixing 0`).
Existing R=0.4 WTA settings, track thresholds and powers 1 and 2 are retained.
The current analyzer uses `L2Relative_AK4PF_OO_V1.txt` for nominal OO JEC;
Run-3 Jet ID and the OO whole-event veto remain enabled. The requested
resources remain 120 minutes and 2,500 MB per job.
These resource settings still require a real-data runtime measurement with mixing.

Rebuild the worker package from the current checkout with:

```bash
python3 makeOOAnalysisTar.py
```

This creates `eecOO.tar.gz` with current sources, calibration payloads, veto
maps and mixing lists. Cards, executables and object files are excluded.
`worker_manifest.json` inside the archive records packaged-file hashes, the
Git revision and whether the source checkout contains local changes. The
submitted card is supplied separately by CRAB.

**Baseline counting policy:** the part01 pilot intentionally retains the
existing ordinary OO/PbPb counting and matching behavior, as requested by
the owner, for comparison with a later counting update. During v1 preparation,
configuration checks, the clean CMSSW build and worker failure-handling
regressions passed. The current configuration enables Vanderbilt, Purdue
and Nebraska.

The ordinary data path uses the legacy branch entry count for the donor chain,
while the `OOClosurePilot 1` MC path uses the total-chain count. The existing
two-file regression demonstrates 2 versus 5 entries for files with 2 and 3
events. The ordinary hiBin matcher retains its existing window widening and
32-bit self-event check. No analyzer changes were made for this baseline.

The v1 task used 14 files/job (950 jobs for part01) and was killed after empty
signal forests caused repeated initialization failures. Its frozen configuration,
card, worker, input lists and archive copies, hashes, submission details and
failure evidence remain under
`submissionRecords/eecAnalysis_OO_mixedPerp_part01_Run3_v1/`.

## ppRef mixed/perpendicular setup: one signal file per job

`crabEECppref5v36.py` uses the same 10,000-job task cap and worker settings
as OO, with ppRef tasks filled to that cap.
The source `ppData2024ForestFiles.txt` contains 31,179 unique forests. The
35 donors in `../mixingFileList/ppRef2024_pool_pilot.txt` are reserved from
that source and excluded from every signal partition. Donors are selected
at evenly spaced midpoint indices across the ordered source, covering its
five PPRefHardProbes datasets in proportion to their file counts.

The remaining 31,144 files retain source-list order and fill four consecutive
tasks: parts 01–03 each contain 10,000 files/jobs, and part04 contains 1,144.
Every part uses `chunk01`; one file per job makes the 10,000-file cap the
maximum configured number of jobs per task. Lists are under
`ppRefForestfiles/ppRefForestFiles_partNN_chunk01.txt`, with `NN` from 01 to 04.
The manifest records source/pool hashes, donor selection and task counts.
The previous ten round-robin lists, their manifest and validation report,
and snapshots of their splitter and CRAB configuration are preserved under
`ppRefForestfiles/previous_10parts/`. The reserved pool, analysis card and
worker payload are unchanged by this repartitioning.
From `crab/`, reproduce or verify the setup and rebuild the shared worker:

```bash
python3 splitPPRefForestFiles.py
python3 makeOOAnalysisTar.py
```

Identical outputs are verified without rewriting them; changed outputs are
rejected. The CRAB configuration also rejects duplicate signal files,
overlap with the pool, and tasks larger than 10,000 one-file jobs.

`cardEECpp5v36TeV.input` selects mixed plus perpendicular cones
(`BackgroundMethods 6`), regular donor forests (`MegaSkimMixing 0`) and the
explicit ppRef pool. It retains unbinned option-0 matching, which acts as
Vz-only matching for ppRef because its reader supplies a constant hiBin.
The card keeps the existing R=0.4 WTA selection, 0.5 GeV track threshold,
EEC binning, powers 1 and 2, ppRef corrections, Jet ID and event veto.
Covariance calculation is disabled as in the OO CRAB card.

**Donor provenance:** these files come from the supplied HardProbes-filtered
forests. They are not an independently established minimum-bias/ZeroBias
pool. `TriggerSelection 8` selects ZeroBias for signal events; the mixing
reader disables its trigger selection for donors. This setup preserves the
existing ordinary OO/ppRef mixer, including its legacy donor-chain entry
count and event-identity comparison described above. File-level separation
does not establish absence of duplicate physical events across forests.

The default task is part01, chunk01. Select parts 1–4 with `PPREF_CRAB_PART`;
`PPREF_CRAB_CHUNK` must be 1. With the usual CMSSW/CRAB environment and CMS
proxy initialized, submit that task from `crab/` with:

```bash
PPREF_CRAB_PART=1 PPREF_CRAB_CHUNK=1 crab submit -c crabEECppref5v36.py
```

To submit all four prepared tasks:

```bash
for part in {1..4}; do
  PPREF_CRAB_PART=$part PPREF_CRAB_CHUNK=1 crab submit -c crabEECppref5v36.py || break
done
```

Names include the partition, for example
`eecAnalysis_ppRef_mixedPerp_part01_chunk01_Run3_v2`. Each job processes one
signal file, requests 120 minutes and 2,500 MB, and allows Vanderbilt, Purdue
or Nebraska. Outputs go to `T2_CH_CERN` under
`/store/group/phys_heavyions/jvelazqu/eec/<requestName>/`.
`eecOO.tar.gz` is the shared OO/ppRef worker archive and includes the ppRef
mixing pool; the separately submitted card selects the collision system.
Preparation does not submit CRAB tasks. The original donor and worker
validation is retained for the unchanged pool and payload: all 35 donors
passed remote metadata checks for nonempty, aligned trees and required
branches; the packaged analysis compiled under CMSSW_15_0_11, and all five
worker regressions passed. That evidence remains in
`ppRefForestfiles/previous_10parts/setup_validation.json` and
`ppRefForestfiles/pool_metadata.json`. The report for the current
partition/configuration checks is `ppRefForestfiles/setup_validation.json`.
These checks do not establish matching correctness, full-job runtime or
physics closure.

## Example commands

1. Create local tar ball:
    ```
    ./makeEECAnalysisTar.sh
    ```

2. Move the tar ball to lxplus:
   ```
   scp eec5TeV.tar.gz jviinika@lxplus.cern.ch:/afs/cern.ch/work/j/jviinika/private/eecAnalysis/CMSSW_12_4_8/src/eecAnalysiss
   ```

3. ssh to lxplus and set up the environment:
   ```
   ssh jviinika@lxplus.cern.ch
   cd /afs/cern.ch/work/j/jviinika/private/eecAnalysis/CMSSW_12_4_8/src/eecAnalysis
   cmsenv
   voms-proxy-init --voms cms
   ```

4. Modify the input card and crab configuration at lxplus:
   ```
   vim cardEECPbPb.input
   vim crabEECPbPb.py
   ```

5. Send the crab jobs:
   ```
   crab submit -c crabEECPbPb.py
   ```

## Executables

**compileAndRun.sh**

This is the script that is executed by CRAB. What is does is that it unzips the analysis tar ball, compiles the code, and runs the analysis.

## Crab configuration files

The configuration files given here are configured to sace the histograms to my (jviinika) EOS space. You will need to change that before you can use these. To do this, just put to variable config.Data.outLFNDirBase a path to a folder to which you have write right and possibly also change the config.Site.storageSite to a site to which to can write. Otherwise for each run you should set the variable infoString, which is the name of the run. It should describe the current configuration used. If you do not know where the input files for the analysis are located, setting fileLocation variable to '2' searches for the files. It is slightly slower to execute, but works with files located everywhere.

**crabEECPbPb.py** /
**crabEECPp.py**

Configuration file for running over the 2018 PbPb or 2017 pp dataset.

**crabEECPbPbMC.py** /
**crabEECPpMC.py**

Configuration file for running over the Pythia+Hydjet or Pythia8 samples.

## Analysis configuration files

These are the operational copies in `crab/`. Their paths are unchanged by
the Run-2 archive move. The Run-3 R=0.4 cards enable HIN Jet ID and the
whole-event veto as described in the [main README](../README.md). The CRAB configurations use local
card basenames, so resolve them from this directory or a staging directory
containing these copies. The archived root-level Run-2 local cards in
[`../cards/archive/run2/`](../cards/archive/run2/README.md) are independent
configurations; do not redirect the CRAB configurations to that archive.
See the [repository card inventory](../cards/README.md).

The configuration for the analysis is provided in card.input files. For example the collision system, data vs. MC running, and several analysis cuts can be defined in this file. The files contain comments on what each variable mean, so it should be straightforward to modify to your needs. The most important example files are.

**cardEECPbPb.input** /
**cardEECPp.input**

Configuration for running with 2018 PbPb or 2017 pp data.

**cardEECPbPbMC.input** /
**cardEECPpMC.input**

Configuration for running over the Pythia+Hydjet or Pythia8 samples.

## File lists

File lists are needed by the CRAB configuration. They contain names for all the files that are analyzed in the configured run. The following file lists are included:

**PbPbData2018_flowJetSkims_80and100triggers.txt**

These are the default files for the 2018 PbPb data. These are forested from /HIHardProbes/HIRun2018A-04Apr2019-v1/AOD primary dataset and either calo jet 80 or calo jet 100 triggers are required for each event in this skim.

**ppData2017_HighEGJet.txt**

Default file list for 2017 pp data. The files are forested from /HighEGJet/Run2017G-17Nov2017-v2/AOD primary dataset.

**newPbPbMcSkim.txt**

Default file list for Pythia+Hydjet simulation. The files are forested from /DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM dataset.

**ppMC2017_mergedFiles.txt**

Default file list for Pythia8 simulation. The files are forested from /QCD_pThat-15_Dijet_TuneCP5_5p02TeV_pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV_v1-v1/AODSIM dataset.

**triggerEfficiencyFiles2018.txt**

Small sample from /HIHardProbes/HIRun2018A-04Apr2019-v1/AOD dataset with no trigger selection and no tracks for trigger efficiency studies.

**PbPbData2018_MinBiasFiles.txt**

Small minimum bias sample for PbPb2018.

## Other files needed for Crab running

**PSet.py** /
**FrameworkJobReport.xml**

Dummy files. Not used, but need to be there to keep CRAB happy.
