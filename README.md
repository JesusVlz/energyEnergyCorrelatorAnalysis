# Energy-energy correlators in OO and pp collisions at 5.36 TeV

This repository contains the framework used to calculate energy-energy correlators (EECs) inside jets. The current development focuses on Run-3 OO collisions and the pp reference at \(\sqrt{s_{\mathrm{NN}}}=5.36\) TeV. The legacy PbPb/pp analysis used for the published 5.02 TeV result is also retained.

> **Current status:** the reconstructed-level O–O and ppRef workflow is available for local tests, CRAB production, QA, and reflected/perpendicular-cone studies. Track-pair corrections, Run-3 JER, detector response, unfolding, and the final corrected O–O/ppRef ratio are still under validation.

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/JesusVlz/energyEnergyCorrelatorAnalysis.git
cd energyEnergyCorrelatorAnalysis
```

### 2. Compile the analysis and plotting code

```bash
make
make -C plotting
```

Run the executable without arguments to print its command-line help:

```bash
./eecAnalysis
```

### 3. Run a local O–O test

```bash
./eecAnalysis \
  /eos/cms/store/group/phys_heavyions/jvelazqu/OO2025_dataset_MinimumBias_filtered/IonPhysics0/OO_MinimumBias_filtered_IonPhysics0/260908_000643/0000/L1_object_data_541.root \
  cardEECOO.input \
  test_OO.root \
  1 \
  true
```

The same command on one line is:

```bash
./eecAnalysis /eos/cms/store/group/phys_heavyions/jvelazqu/Forest/OO/L1_object_data.root cardEECOO.input test_OO.root 1 true
```

The arguments are:

```text
./eecAnalysis INPUT CARD OUTPUT FILE_LOCATION RUN_LOCAL [MIXING_LIST_INDEX]
```

| Argument | Meaning |
|---|---|
| `INPUT` | One ROOT file or a text file containing one input file per line. In CRAB mode, this is the CRAB job ID. |
| `CARD` | Analysis configuration card. |
| `OUTPUT` | Output ROOT file. |
| `FILE_LOCATION` | `0` Purdue, `1` CERN EOS, `2` Vanderbilt, `3` global xrootd, or `-1` ACCRE/Slurm. For a complete ROOT path or URL, the value is still required but is not prepended. |
| `RUN_LOCAL` | `true` for a local test; `false` for CRAB mode. |
| `MIXING_LIST_INDEX` | Optional mixed-event list index. Normally omitted for local tests. |

### 4. Run a local ppRef test

```bash
./eecAnalysis \
  root://xrootd-vanderbilt.sites.opensciencegrid.org//store/user/jvelazqu/HiForestpp2024/PPRefHardProbes0/ppRef2024_HardProbes_filtered_PPRefHardProbes0/260925_144732/0000/HiForestMiniAOD_1.root \
  cardEECpp5v36TeV.input \
  test_ppRef.root \
  2 \
  true
```

For MC, use `cardEECOOMC.input` or `cardEECpp5v36TeVMC.input` with an O–O or ppRef MC forest.

### 5. Check the output

Open the file and print the configuration saved in it:

```bash
root -l test_OO.root
root -l -b -q 'plotting/checkCard.C("test_OO.root")'
```

Before starting a large production, confirm that:

- the output contains nonzero event, jet, track, and EEC histograms;
- the expected `DataType`, jet radius, WTA axis, trigger, and background method are stored in the card;
- no required tree or branch is missing;
- the JEC and tracking-correction files are found;
- all event, jet, and track weights are finite.

## Run-3 cards

| Sample | `DataType` | Local card | CRAB configuration |
|---|---:|---|---|
| O–O data | 9 | `cardEECOO.input` | `crab/crabEECOO.py` |
| O–O MC | 10 | `cardEECOOMC.input` | `crab/crabEECOOMC.py` |
| ppRef data | 11 | `cardEECpp5v36TeV.input` | `crab/crabEECppref5v36.py` |
| ppRef MC | 12 | `cardEECpp5v36TeVMC.input` | `crab/crabEECpprefMC5v36.py` |

The current baseline uses AK4 jets, \(R=0.4\), the winner-take-all axis, \(|\eta_{\mathrm{jet}}|<1.6\), and \(|\eta_{\mathrm{track}}|<2.4\). O–O uses the 0–30%, 30–50%, and 50–90% centrality intervals, while ppRef is inclusive.

The root-level data cards use `BackgroundMethods 0` for simple local tests. The CRAB data cards use `BackgroundMethods 3` to fill reflected and perpendicular cones. MC baseline cards use `BackgroundMethods 0`.

## Running with CRAB

### 1. Set up CRAB on lxplus

Enter a CMSSW area and initialize the environment:

```bash
cd /path/to/CMSSW_X_Y_Z/src
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --voms cms --valid 192:00
voms-proxy-info --timeleft
```

Then return to the analysis repository.

### 2. Build the analysis tarball

Build a fresh tarball from the exact commit that will be submitted:

```bash
git status --short
git rev-parse HEAD
./makeEECAnalysisTar.sh
cp -f eec5TeV.tar.gz crab/eecOO.tar.gz
```

The copy command is currently necessary because `makeEECAnalysisTar.sh` creates `eec5TeV.tar.gz`, while the Run-3 CRAB configurations and `crab/compileAndRun.sh` expect `eecOO.tar.gz`.

Do not reuse an old tarball. CRAB compiles and runs the source stored inside this tarball, not necessarily the source currently visible in the working directory.

### 3. Review the card and CRAB configuration

```bash
cd crab
vim cardEECOO.input
vim crabEECOO.py
```

For the selected sample, check:

- `jobTag` is new and describes the production;
- `card` and `inputList` select the intended sample;
- `unitsPerJob` is appropriate;
- `outLFNDirBase` points to a writable EOS directory;
- `storageSite` and any site whitelist are correct;
- the card contains the intended data type, trigger, cuts, and background method.

The current input lists are:

| Sample | Input list | Initial files/job |
|---|---|---:|
| O–O data | `OOForestFiles.txt` | 14 |
| O–O MC | `OOMC2025ForestFiles.txt` | 2 |
| ppRef data | `ppData2024ForestFiles.txt` | 5 |
| ppRef MC | `ppMC2024ForestFiles.txt` | 2 |

Use `TriggerSelection 0` for O–O and ppRef MC; no trigger selection is required. For the first MC smoke test, also use `BackgroundMethods 0` and `DisableTrackPairEfficiencyCorrection 1`.

### 4. Submit one pilot job

Before the full submission, temporarily set:

```python
config.Data.totalUnits = config.Data.unitsPerJob
```

Submit the selected configuration:

```bash
crab submit -c crabEECOO.py
```

The other Run-3 submissions are:

```bash
crab submit -c crabEECOOMC.py
crab submit -c crabEECppref5v36.py
crab submit -c crabEECpprefMC5v36.py
```

Validate the pilot output before restoring:

```python
config.Data.totalUnits = len(config.Data.userInputFiles)
```

Use a new `jobTag` for the full production.

### 5. Monitor and control the task

The CRAB project directory normally has the form `JOB_TAG/crab_JOB_TAG`:

```bash
PROJECT=JOB_TAG/crab_JOB_TAG

crab status -d "$PROJECT"
crab status --long -d "$PROJECT"
crab report -d "$PROJECT"
crab getlog --short -d "$PROJECT" --jobids 1
crab getoutput --dump -d "$PROJECT"
```

Resubmit failed jobs:

```bash
crab resubmit -d "$PROJECT"
```

If the job report shows that more runtime or memory is required:

```bash
crab resubmit -d "$PROJECT" --maxjobruntime 180 --maxmemory 3000
```

Stop a task:

```bash
crab kill -d "$PROJECT"
```

## Merge and post-process the output

Merge the successful ROOT files. A simple local example is:

```bash
hadd -f merged_OO.root /path/to/OOfiles/*.root
hadd -f merged_ppRef.root /path/to/ppRefFiles/*.root
```

For large productions, merge in several stages to avoid opening thousands of remote files in one `hadd` command.

Project the sparse histograms:

```bash
./projectHistogramsInSteps.sh merged_OO.root merged_OO_processed.root -w 1
./projectHistogramsInSteps.sh merged_ppRef.root merged_ppRef_processed.root -w 1
```

For files filled with `BackgroundMethods 3`, create explicit background-processing variants:

```bash
# Perpendicular-cone subtraction
./processEnergyEnergyCorrelators.sh merged_OO_processed.root \
  -o merged_OO_perpendicular.root -b 1

# Reflected-cone subtraction
./processEnergyEnergyCorrelators.sh merged_OO_processed.root \
  -o merged_OO_reflected.root -b 2

# Signal without background subtraction
./processEnergyEnergyCorrelators.sh merged_OO_processed.root \
  -o merged_OO_noSubtraction.root -b 3
```

Repeat the same commands for ppRef. Do not use the default `-b 0` unless the file contains a validated mixed-cone background.

Before calculating the O–O/ppRef ratio, verify that both files use the same jet radius, WTA axis, jet and track bins, energy-weight exponent, normalization interval, and processing level.

## Run-3 response and unfolding pilot

`UnfoldingForestReader` supports OO MC (`DataType 10`) and ppRef MC
(`DataType 12`), including reconstructed/generated jet matching and the MiniAOD
track vectors. Required trees, branch types, vector sizes, and event counts are
checked. These paths pass synthetic ROOT regression and response-filling checks;
a production MC forest still needs a one-file validation.

For an initial response pilot, copy the corresponding MC card and set:

```text
McCorrelationType 3
MatchJets 1
TriggerSelection 0
FilledHistograms 223
BackgroundMethods 0
JetUncertainty 0
JetPtWeight 0
SmearDeltaR 0
SmearEnergyWeight 0
DisableTrackPairEfficiencyCorrection 1
```

Mask `223` fills the regular histograms, jet closure, and unfolding response.
Use `479` to include track–particle matching as well. Keep the nominal jet,
track, and centrality bins consistent with the target data. The response uses
generator particles inside reconstructed/generated jet cones regardless of the
accompanying MC EEC correlation mode; it is a jet response for particle EECs,
not a complete correction for single-track and track-pair detector effects.

Run on one MC forest and project the response, for example:

```bash
./eecAnalysis /path/to/OO_MC.root cardOOResponse.input OO_response.root 1 true
./projectHistogramsInSteps.sh OO_response.root OO_response_processed.root -c -r -a -w 1
# Add -m when the pilot card includes track–particle matching.
```

OO/ppRef unfolding requires an explicit processed response file and a positive
iteration count selected through independent closure studies. It cannot select
the legacy pp2017 response automatically. In a ROOT session with RooUnfold
available, the extended macro call is:

```cpp
.L plotting/unfoldEEChistograms.C
unfoldEEChistograms("OO_data_processed.root", "OO_unfolded.root",
                    0, 0, 0, 1, "OO_response_processed.root", chosenIterations);
```

The file and count are used explicitly for the requested split/systematic;
the caller must supply the corresponding studied response. Compatibility checks
cover collision system, MC metadata, energy exponent, jet type/axis/radius,
binning, and response content. All OO centralities are retained. For Run-3, the
macro selects complete native card bins inside its 120–200 GeV study window
(120–160 and 160–200 GeV for the current cards).

No dedicated Run-3 JER scale factors or EEC smearing files are installed.
Nonzero Run-3 `JetUncertainty`, `JetPtWeight`, or smearing switches fail before
initialization. Nominal pilot matrices contain the MC detector resolution;
the additional data/MC resolution calibration, independent closure, iteration
study, and a RooUnfold end-to-end check remain necessary. Existing processed
OO/ppRef MC files in this workspace do not contain response/jet-closure outputs.

## Current progress

| Component | Status |
|---|---|
| O–O and ppRef data types 9–12 | Implemented in the main analysis path |
| AK4, \(R=0.4\), WTA cards | Configured for O–O and ppRef data/MC |
| O–O and ppRef JEC | Selected in the analyzer; closure is still required |
| O–O tracking correction | Implemented with the 2025 nominal table |
| ppRef tracking correction | Nominal table included with checksum and loader checks |
| RecoReco production and QA | Current baseline workflow |
| RecoGen and GenGen | Routing added; branch-by-branch validation is still required |
| Reflected and perpendicular cones | Available with `BackgroundMethods 3` |
| Mixed events | O–O pilot only; not part of the nominal baseline |
| Track-pair efficiency | Disabled/unity until a Run-3 correction is derived |
| Run-3 JER | Dedicated inputs outstanding; unsupported variations fail clearly |
| Response and unfolding | Run-3 reader and explicit response configuration implemented; synthetic checks pass; production closure outstanding |
| Final corrected O–O/ppRef ratio | Not yet final |

## Next validation steps

Do these approximately in this order:

1. Rebuild the CRAB tarball and run one-file RecoReco tests for all four samples.
2. Validate tracking closure; OO keeps the common PbPb EEC coordinate convention.
3. Validate the installed O–O/ppRef marginal event weights in an event-level rerun and establish generator normalization.
4. Validate RecoGen and GenGen trees and branches.
5. Perform jet matching, JEC closure, and Run-3 JER studies.
6. Perform single-track efficiency and fake-rate closure.
7. Derive and validate the residual track-pair efficiency.
8. Validate reflected/perpendicular background subtraction and its normalization.
9. Produce dedicated O–O/ppRef response files, select iterations with independent closure, and validate RooUnfold end to end.
10. Build the Run-3 systematic uncertainties and then produce the corrected O–O/ppRef ratio.

## Repository structure

| Path | Contents |
|---|---|
| `eecAnalysis.cxx` | Main executable and command-line interface |
| `src/` | Forest readers, analyzer, histograms, and correction code |
| `plotting/` | Projection, processing, QA, comparison, and unfolding macros |
| `crab/` | CRAB configurations, cards, input lists, and wrapper scripts |
| `jetEnergyCorrections/` | Jet-energy correction files |
| `trackCorrectionTables/` | Single-track and track-pair correction tables |
| `mixingFileList/` | Mixed-event input lists |
| `smearingFiles/` | Detector-response smearing inputs |
| `skimmer/` | Minimum-bias skimming tools for event mixing |
| `toySimulation/` | Uniform-particle toy simulations |

## Legacy PbPb/pp paper plots

The legacy part of this repository reproduces the figures from *Observation of nuclear modification of energy-energy correlators inside jets in heavy ion collisions*, [Phys. Lett. B 866 (2025) 139556](https://doi.org/10.1016/j.physletb.2025.139556), [arXiv:2503.19993](https://arxiv.org/abs/2503.19993).

On lxplus, with Vanderbilt T2 read access:

```bash
voms-proxy-init --voms cms
xrdcp root://xrootd-vanderbilt.sites.opensciencegrid.org//store/user/jviinika/energyEnergyCorrelatorAnalysisFiles/eecAllFinalResultFiles.tar.gz .
tar -xvzf eecAllFinalResultFiles.tar.gz
./makeAllPaperPlots.sh
```

The publication figures are written under `paperFigures/`.

