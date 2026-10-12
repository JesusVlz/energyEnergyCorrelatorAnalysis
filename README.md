# Energy-energy correlators in O–O and ppRef collisions at 5.36 TeV

This repository contains the analysis framework used to calculate energy-energy correlators (EECs) inside jets. The current development focuses on Run-3 oxygen–oxygen (O–O) collisions and the 5.36 TeV proton–proton reference (ppRef), while retaining the legacy PbPb/pp analysis used for the published 5.02 TeV measurement.

> **Development status:** the Run-3 detector-level workflow is under active validation. Reconstructed-level data/MC production, QA, and reflected/perpendicular-cone studies are the present focus. Track-pair corrections, Run-3 JER, response production, unfolding, and the final corrected O–O/ppRef ratio are not yet production-ready.
>
> This status summary corresponds to repository commit `f57c623` (October 9, 2026). Run `git log -1 --oneline` and review newer changes before launching a large production.

## Analysis scope and baseline configuration

The Run-3 data types are:

| `DataType` | Collision system | Sample type | Main card |
|---:|---|---|---|
| 9 | O–O 2025 | Data | `cardEECOO.input` |
| 10 | O–O 2025 | MC | `cardEECOOMC.input` |
| 11 | ppRef 5.36 TeV | Data | `cardEECpp5v36TeV.input` |
| 12 | ppRef 5.36 TeV | MC | `cardEECpp5v36TeVMC.input` |

The nominal Run-3 comparison currently uses:

- AK4 jets with `JetRadius 0.4`;
- winner-take-all axes with `JetAxis 1`;
- jet acceptance `|eta_jet| < 1.6`;
- track acceptance `|eta_track| < 2.4`;
- reconstructed tracks down to 0.5 GeV in the analysis files;
- the main first-look EEC selection `pT_ch > 2 GeV`, obtained during projection/plotting;
- energy-weight exponents `n = 1` and `n = 2`;
- O–O centrality intervals 0–30%, 30–50%, and 50–90%;
- inclusive ppRef, represented by the ppRef `C0` bin;
- `HLT_MinimumBiasHF_OR_BptxAND_v1` for O–O data;
- `HLT_PPRefZeroBias_v6` for ppRef data.

For detector-level background studies, `BackgroundMethods 3` fills both reflected and perpendicular cones. Mixed events are not part of the nominal Run-3 baseline.

## Repository layout

| Path | Purpose |
|---|---|
| `eecAnalysis.cxx` | Main executable and command-line interface |
| `src/` | Event readers, analyzer, histogram definitions, corrections, and response code |
| `plotting/` | Projection, post-processing, QA, comparison, closure, and unfolding macros |
| `crab/` | CRAB configurations, Run-3 cards, forest file lists, and CRAB wrapper scripts |
| `jetEnergyCorrections/` | Jet-energy correction files |
| `trackCorrectionTables/` | Single-track and track-pair correction tables |
| `mixingFileList/` | Mixed-event input lists, including the current O–O pilot list |
| `smearingFiles/` | Legacy detector-response smearing inputs |
| `skimmer/` | Minimum-bias skimming tools used for event mixing |
| `toySimulation/` | Uniform-particle toy simulations |

The root-level cards are intended for local tests and controlled analysis jobs. CRAB transfers the corresponding copies from `crab/`. Always compare the two copies before production because they can diverge.

## Requirements

For local compilation and analysis:

- ROOT with `root-config` available;
- a C++ compiler supported by the installed ROOT version;
- GNU Make;
- network or local access to the input HiForest files.

For CRAB production on lxplus:

- a working CMSSW area;
- the CRAB client environment;
- a valid CMS VOMS proxy;
- write access to the configured output storage site and EOS directory.

## Clone and compile

```bash
git clone https://github.com/JesusVlz/energyEnergyCorrelatorAnalysis.git
cd energyEnergyCorrelatorAnalysis

git log -1 --oneline
make
make -C plotting
```

Run the executable without arguments to print its command-line help:

```bash
./eecAnalysis
```

## Run a local test

The executable has the following interface:

```text
./eecAnalysis INPUT CARD OUTPUT FILE_LOCATION RUN_LOCAL [MIXING_LIST_INDEX]
```

| Argument | Meaning |
|---|---|
| `INPUT` | One ROOT file or a text file containing one input path per line. In CRAB mode this is the CRAB job ID. |
| `CARD` | Analysis configuration card. |
| `OUTPUT` | Output ROOT histogram file. |
| `FILE_LOCATION` | `0` Purdue, `1` CERN EOS, `2` Vanderbilt, `3` global xrootd, or `-1` ACCRE/Slurm mode. |
| `RUN_LOCAL` | `true` for a direct ROOT file or ordinary text list; `false` for CRAB mode. |
| `MIXING_LIST_INDEX` | Optional mixed-event list index. The default is `-1`, which derives it from the CRAB job ID. |

When `RUN_LOCAL=true` and `INPUT` is already a complete local path, `/eos/cms/...` path, or `root://...` URL, the program uses that path directly. `FILE_LOCATION` is still required syntactically but is not prepended to a direct ROOT filename.

### O–O data example

```bash
./eecAnalysis \
  /eos/cms/store/group/phys_heavyions/jvelazqu/OO2025_dataset_MinimumBias_filtered/IonPhysics0/OO_MinimumBias_filtered_IonPhysics0/260908_000643/0000/L1_object_data_541.root \
  cardEECOO.input \
  test_OO.root \
  1 \
  true
```

This is equivalent to the one-line command:

```bash
./eecAnalysis /eos/cms/store/group/phys_heavyions/jvelazqu/OO2025_dataset_MinimumBias_filtered/IonPhysics0/OO_MinimumBias_filtered_IonPhysics0/260908_000643/0000/L1_object_data_541.root cardEECOO.input test_OO.root 1 true
```

### Run from a text file list

Create a list containing one complete path or URL per line:

```text
/eos/cms/store/group/.../forest_1.root
/eos/cms/store/group/.../forest_2.root
```

Then run:

```bash
./eecAnalysis testFileList.txt cardEECOO.input test_OO_list.root 1 true
```

For a ppRef test, use `cardEECpp5v36TeV.input`. For MC, use `cardEECOOMC.input` or `cardEECpp5v36TeVMC.input`.

### Inspect the configuration saved in the output

Every output contains an `EECCard`. Print it with:

```bash
root -l -b -q 'plotting/checkCard.C("test_OO.root")'
```

Before accepting a test, verify:

- the expected `DataType`, `JetType`, `JetAxis`, and `JetRadius`;
- the expected trigger and background bit masks;
- nonzero accepted event, jet, track, and EEC counts;
- no missing-tree or missing-branch warnings;
- finite JEC and tracking weights;
- the Git hash printed/stored by the production tarball;
- both perpendicular-cone orientations when `BackgroundMethods 3` is used.

## Important card parameters

### MC correlation modes

| `McCorrelationType` | Jet level | Constituent level | Intended use |
|---:|---|---|---|
| 0 | Reconstructed | Reconstructed | Nominal detector-level distribution |
| 1 | Reconstructed | Generator | Reco/gen constituent comparison |
| 2 | Generator | Reconstructed | Specialized cross-check |
| 3 | Generator | Generator | Particle-level reference |

### Histogram bit mask

`FilledHistograms` is the sum of the desired bits:

| Bit value | Histograms |
|---:|---|
| 1 | Event information |
| 2 | Jets |
| 4 | Tracks |
| 8 | Multiplicity and density inside jets |
| 16 | Nominal EECs |
| 32 | EEC systematic variations |
| 64 | Jet-pT closure |
| 128 | Jet-pT unfolding response |
| 256 | Track/particle matching |

The regular detector-level configuration is `31`. Do not automatically use `511`: dedicated cards are easier to validate and require less memory.

### Background bit mask during filling

| Bit value | Method |
|---:|---|
| 1 | Reflected-eta cone |
| 2 | Perpendicular cones at `+pi/2` and `-pi/2` |
| 4 | Mixed events |

Therefore:

- `BackgroundMethods 0`: no cone background; use for initial MC closure;
- `BackgroundMethods 3`: reflected plus perpendicular cones; use for Run-3 data background studies;
- do not enable bit 4 until the mixed-event sample and matching are validated.

The post-processing background indices are different from the card bit mask:

| `processEnergyEnergyCorrelators.sh -b` | Post-processing choice |
|---:|---|
| 0 | Mixed-cone subtraction |
| 1 | Perpendicular-cone subtraction |
| 2 | Reflected-cone subtraction |
| 3 | No subtraction |
| 4 | MC truth-level signal/background separation |

Do not use the script default `-b 0` for a Run-3 file that does not contain validated mixed cones.

## CRAB production

The current Run-3 CRAB configurations are:

| Sample | CRAB configuration | Input list | Current files/job |
|---|---|---|---:|
| O–O data | `crabEECOO.py` | `OOForestFiles.txt` | 14 |
| O–O MC | `crabEECOOMC.py` | `OOMC2025ForestFiles.txt` | 2 |
| ppRef data | `crabEECppref5v36.py` | `ppData2024ForestFiles.txt` | 5 |
| ppRef MC | `crabEECpprefMC5v36.py` | `ppMC2024ForestFiles.txt` | 2 |

These numbers are starting points, not universal recommendations. Re-evaluate them after one pilot submission using the observed runtime, memory, and remote-read performance.

### 1. Set up the lxplus environment

From a CMSSW area:

```bash
cd /path/to/CMSSW_X_Y_Z/src
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --voms cms --valid 192:00
voms-proxy-info --timeleft
```

Then enter the analysis repository.

### 2. Rebuild the analysis tarball from the exact commit

Do not submit the old committed `crab/eecOO.tar.gz`. Build a clean tarball from the checkout you intend to run:

```bash
git status --short
git rev-parse HEAD
./makeEECAnalysisTar.sh
cp -f eec5TeV.tar.gz crab/eecOO.tar.gz
```

This rename is currently required because `makeEECAnalysisTar.sh` creates `eec5TeV.tar.gz`, while the four Run-3 CRAB configurations request `eecOO.tar.gz`.

Confirm that the tarball does not contain cards, compiled object files, or an old executable:

```bash
tar -tzf crab/eecOO.tar.gz | grep -E '(^|/)card.*\.input$|\.o$|(^|/)eecAnalysis$'
```

The expected result is no output. Cards are transferred separately by CRAB, which prevents a card inside the tarball from silently overwriting the intended production card.

### 3. Review the production card and CRAB configuration

Move to the CRAB directory:

```bash
cd crab
```

For the chosen sample, check:

- a new, descriptive `jobTag`/`requestName`;
- the intended `card` and `inputList`;
- `unitsPerJob`;
- `fileLocation` and the site whitelist;
- `outLFNDirBase` and `storageSite`;
- the card's `DataType`, trigger, jet axis/radius, cuts, and background mask;
- that the root-level and CRAB copies of the card agree where they should.

Example card comparison:

```bash
diff -u ../cardEECOO.input cardEECOO.input
```

Some differences may be intentional—for example, a root-level no-background smoke-test card versus a CRAB data-cone card—but every difference should be understood before submission.

For the first MC validation job, use a dedicated copy of the MC card with:

```text
TriggerSelection 0
BackgroundMethods 0
DisableTrackPairEfficiencyCorrection 1
JetUncertainty 0
JetPtWeight 0
SmearDeltaR 0
SmearEnergyWeight 0
SkipCovarianceMatrix 1
```

Only restore `TriggerSelection 8` after verifying that the relevant MC trigger branch exists and has meaningful acceptance.

### 4. Perform a one-job dry run

Before submitting thousands of jobs, temporarily set:

```python
config.Data.totalUnits = config.Data.unitsPerJob
```

Submit the appropriate configuration:

```bash
crab submit -c crabEECOO.py
```

Other choices are:

```bash
crab submit -c crabEECOOMC.py
crab submit -c crabEECppref5v36.py
crab submit -c crabEECpprefMC5v36.py
```

CRAB prints the project directory after submission. With the current configuration layout it is normally:

```text
JOB_TAG/crab_JOB_TAG
```

### 5. Monitor and debug the task

Set the project path once:

```bash
PROJECT=JOB_TAG/crab_JOB_TAG
```

Then use:

```bash
crab status -d "$PROJECT"
crab status --long -d "$PROJECT"
crab report -d "$PROJECT"
crab getlog --short -d "$PROJECT" --jobids 1
crab getoutput --dump -d "$PROJECT"
```

For failed jobs:

```bash
crab resubmit -d "$PROJECT"
```

If the failure is caused by insufficient runtime or memory, first inspect the job report and then resubmit with justified limits, for example:

```bash
crab resubmit -d "$PROJECT" --maxjobruntime 180 --maxmemory 3000
```

To stop a task:

```bash
crab kill -d "$PROJECT"
```

To monitor several tasks, place one `jobTag` per line in a text file and use:

```bash
./checkCrabJobs.sh submittedJobs.txt
./checkCrabJobs.sh submittedJobs.txt -r
./checkCrabJobs.sh submittedJobs.txt -r -t 180
```

The helper assumes the project layout `JOB_TAG/crab_JOB_TAG`.

### 6. Validate the pilot output before full submission

For each sample, inspect the stdout/stderr and output ROOT file. Confirm:

- the correct Git hash and card are printed;
- the correction files are found;
- all required trees and branches are connected;
- event weights are finite and the event count is nonzero;
- expected trigger and event selections are applied;
- the jet collection is AK4 and the axis is WTA;
- the output can be opened and projected;
- data cone jobs contain both reflected and perpendicular pairings.

Only then restore `config.Data.totalUnits = len(config.Data.userInputFiles)` and submit the full task with a new `jobTag`.

## Merge, project, and process the output

CRAB writes one ROOT file per successful job to the configured EOS destination. For small tests, obtain the output PFNs with:

```bash
crab getoutput --dump -d "$PROJECT"
```

For a large O–O production, merge in stages rather than giving thousands of remote files to one `hadd` process. The `crab/addHistogramsInSteps.sh` helper performs staged merging, but it currently hard-codes the FNAL EOS endpoint; check or update the endpoint before using it for files stored at CERN.

After merging, project the sparse histograms:

```bash
./projectHistogramsInSteps.sh merged_OO.root merged_OO_processed.root -w 1
./projectHistogramsInSteps.sh merged_ppRef.root merged_ppRef_processed.root -w 1
```

Print the saved cards:

```bash
root -l -b -q 'plotting/checkCard.C("merged_OO_processed.root")'
root -l -b -q 'plotting/checkCard.C("merged_ppRef_processed.root")'
```

For a file filled with `BackgroundMethods 3`, produce explicit post-processing variants:

```bash
./processEnergyEnergyCorrelators.sh merged_OO_processed.root \
  -o merged_OO_perpendicular.root -b 1

./processEnergyEnergyCorrelators.sh merged_OO_processed.root \
  -o merged_OO_reflected.root -b 2

./processEnergyEnergyCorrelators.sh merged_OO_processed.root \
  -o merged_OO_noSubtraction.root -b 3
```

Repeat the same commands for ppRef. Do not form the O–O/ppRef ratio until both files use identical jet/track bins, radius, axis, weight exponent, normalization interval, and processing level.

Additional projection modes include:

```bash
# Jet closure
./projectHistogramsInSteps.sh jetClosure.root jetClosure_processed.root -n -c

# Jet/EEC response
./projectHistogramsInSteps.sh response.root response_processed.root -n -r -a

# Track/particle matching
./projectHistogramsInSteps.sh matching.root matching_processed.root -n -m

# Tracking and pair-correction systematic histograms
./projectHistogramsInSteps.sh variations.root variations_processed.root -n -e
```

Use these only after the corresponding Run-3 reader and correction path has been validated.

## Current Run-3 implementation status

| Component | Current status |
|---|---|
| Run-3 data types 9–12 | Implemented in the main `ForestReader`/`HighForestReader` path |
| AK4 O–O and ppRef cards | Aligned to `R=0.4` and WTA for the baseline |
| O–O JEC | Configured with `L2Relative_AK4PF_OO_V1.txt` |
| ppRef JEC | Run-3 relative/residual files are selected explicitly |
| O–O tracking | Nominal 2025 efficiency table is present; the lookup is capped at 6.1 GeV |
| ppRef tracking | Analyzer class exists, but the referenced ROOT table is absent from the repository and must be deployed reproducibly |
| Generator-level reader | Explicit O–O MC and ppRef MC routing has been added; one-file RecoGen/GenGen branch validation is still required |
| Reflected/perpendicular cones | Available through `BackgroundMethods 3`; used in dedicated data cards |
| O–O mixed events | Pilot list and selection logic exist; not part of the nominal validated workflow |
| Track-pair efficiency | Forced off/unity for O–O and ppRef until dedicated residual corrections are derived |
| Run-3 JER | No dedicated O–O/ppRef JER scale-factor configuration yet |
| Unfolding reader | Still supports only the legacy data types 0–3; Run-3 response and matching jobs remain blocked |
| Run-3 unfolding configuration | Dedicated O–O/ppRef response files and closure are not implemented |
| Final O–O/ppRef ratio | Not yet final; background, response, closure, and systematic validation remain |

### Current MC-weight behavior

The present code should be treated as a diagnostic configuration:

- O–O MC applies the currently installed vertex-z and centrality polynomials;
- O–O MC currently forces the forest event/pThat weight to unity;
- ppRef MC currently uses unity vertex-z, centrality, and forest event weights in the event loop;
- `JetPtWeight` must remain zero until an optional truth/reco jet-spectrum weight is derived and closed.

Validate each factor independently before combining weighted MC samples.

## Prioritized validation and development checklist

### P0 — Quick blockers before new production

- [ ] Rebuild `crab/eecOO.tar.gz` from the exact intended Git commit for every production.
- [ ] Resolve the `eec5TeV.tar.gz` versus `eecOO.tar.gz` naming mismatch permanently.
- [ ] Install and version the ppRef 2024 tracking-efficiency ROOT table referenced by `TrkEff2024ppRef`.
- [ ] Use detector eta, not jet-relative eta, when evaluating single-track corrections inside EEC pair construction.
- [ ] Add the missing Condor/CRAB ignore rules and untrack generated logs, caches, and task directories.
- [ ] Synchronize or clearly document intentional differences between root-level and `crab/` card copies.
- [ ] Correct the hard-coded `5.02 TeV` label in `plotting/EECHistogramManager.cxx` for O–O and ppRef.
- [ ] Run one-file RecoReco smoke tests for O–O data, O–O MC, ppRef data, and ppRef MC.

Recommended ignore rules:

```gitignore
/condor/output_condor/
/condor/error_condor/
/condor/log_condor/
/crab/*/crab_*/
/crab/__pycache__/
*.pyc
/crab.log
/crab/crab.log
```

### P1 — MC readers, weights, and detector closure

- [ ] Validate every O–O/ppRef generator tree and branch with `McCorrelationType 1` and `3`.
- [ ] Extend `UnfoldingForestReader` to data types 10 and 12 without fallback to legacy pp.
- [ ] Confirm the Run-3 raw-jet branch convention (`rawpt` versus `jtptUncorrected`) for each current forest.
- [ ] Make Run-3 MC event weights configurable; validate pThat, vertex-z, and O–O centrality weights one at a time.
- [ ] Run reconstructed-to-generator jet matching and JEC closure.
- [ ] Implement and validate O–O/ppRef JER corrections and variations.
- [ ] Run single-track efficiency/fake/secondary/duplicate closure versus track pT, eta, jet pT, and O–O centrality.
- [ ] Validate the O–O minimum-bias and ppRef zero-bias triggers with an unbiased denominator.

Suggested small MC job matrix:

| Purpose | `McCorrelationType` | `MatchJets` | `FilledHistograms` |
|---|---:|---:|---:|
| Detector-level baseline | 0 | 0 | 31 |
| Reco/gen constituent comparison | 1 | 0 | 31 |
| Particle-level reference | 3 | 0 | 31 |
| Matched-jet study | 1 | 1 | 31 |
| Jet closure | 0 | 1 | 65 |
| Unfolding response | 0 | 0 | 129 |
| Track/particle matching | 0 | 0 | 257 |

The final two modes remain blocked until `UnfoldingForestReader` is Run-3 capable.

### P2 — Pair and background validation

- [ ] Derive residual track-pair efficiency after applying single-track weights.
- [ ] Produce independent O–O and ppRef pair-correction files with the histogram names expected by `TrackPairEfficiencyCorrector`.
- [ ] Validate reflected and averaged perpendicular cones versus centrality, jet pT, track pT, and delta-R.
- [ ] Derive Run-3 background normalization scales instead of falling back to Run-2 pp/PbPb values.
- [ ] Perform MC background-subtraction closure.
- [ ] Use mixed events only if cone methods fail closure or leave a significant residual.

### P3 — Response, unfolding, and final result

- [ ] Add explicit O–O centrality-dependent and ppRef-inclusive unfolding configurations.
- [ ] Construct response matrices from one MC split and close on an independent split.
- [ ] Perform self-closure, cross-closure, refolding, iteration scans, and prior-dependence tests.
- [ ] Derive MC non-closure separately for O–O and ppRef.
- [ ] Add only genuine Run-3 systematic variations: JEC/JER, tracking, pair efficiency, background, trigger, unfolding, normalization, and non-closure.
- [ ] Unfold O–O and ppRef independently, subtract the validated background, normalize over the same delta-R interval, and then form each O–O centrality/ppRef `C0` ratio.

## Expected final analysis sequence

1. Apply JEC and validated single-track corrections during histogram filling.
2. Apply a validated residual pair correction during filling, if required.
3. Project same-jet and background-cone histograms.
4. Unfold O–O and ppRef independently.
5. Construct and subtract the validated background after unfolding.
6. Normalize both systems over the same delta-R interval.
7. Form O–O centrality divided by inclusive ppRef.
8. Propagate statistical uncertainties and correlated/uncorrelated systematics.

## Legacy PbPb/pp paper plots

The original framework can reproduce the figures from:

> *Observation of nuclear modification of energy-energy correlators inside jets in heavy ion collisions*, Phys. Lett. B 866 (2025) 139556, [doi:10.1016/j.physletb.2025.139556](https://doi.org/10.1016/j.physletb.2025.139556), [arXiv:2503.19993](https://arxiv.org/abs/2503.19993).

On lxplus, with Vanderbilt T2 read access:

```bash
voms-proxy-init --voms cms
xrdcp root://xrootd-vanderbilt.sites.opensciencegrid.org//store/user/jviinika/energyEnergyCorrelatorAnalysisFiles/eecAllFinalResultFiles.tar.gz .
tar -xvzf eecAllFinalResultFiles.tar.gz
./makeAllPaperPlots.sh
```

The reproduced publication figures are written under `paperFigures/`.

## Contributing

Before committing:

```bash
git status --short
git diff --check
git diff
```

Do not commit CRAB task directories, Condor stdout/stderr/log files, Python bytecode, merged ROOT outputs, compiled objects, or production tarballs. Keep each physics change accompanied by a small reproducible test and record the input forest, card, Git hash, and expected output checks.

