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

### Optional QA histogram binning

Add the following setting to an analysis card when the inclusive QA spectra
should use the same pT boundaries as the EEC selections:

```text
MatchEECBinning 1
```

This uses `JetPtBinEdgesEEC` for inclusive-jet pT and reflected-cone jet pT,
and `TrackPtBinEdgesEEC` for corrected and uncorrected track pT. It makes
comparisons with the selected EEC bins direct. Omitting the setting, or setting
it to `0`, keeps the established QA binning; `ConfigurationCard::HasKey` allows
older cards to omit this optional key. The maximum-particle-pT and jet-closure
axes retain their existing binning because their consumers apply separate
particle thresholds and jet-pT selections that require those axes.

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

## Plotting and OO diagnostics

Build `plotting/DrawingClasses.so` with `make -C plotting`, then run the macros
from the repository root in a ROOT environment. Pass the input paths explicitly;
the filenames in the examples below are placeholders for your projected or raw
files. Outputs normally go to `figures/`.

The general `plotEEChistograms.C` configuration now selects raw same-jet EECs
instead of track QA, with centrality bin `0`, jet bin `0`, track threshold bin
`3`, and the subevent drawing style enabled for OO MC. In the current OO cards,
these indices mean 0–30%, 40–80 GeV jets, and track pT above 2 GeV. Its default
figure suffix is `_OO_C0-30_GenGen`; edit the selection switches and suffix when
using another sample. OO MC legends identify Pythia+Hijing at 5.36 TeV and use
Hijing labels for the background subevents. Subevent figures use one legend and
choose their vertical range from all displayed curves and statistical errors,
so small components are visible without a fixed legacy y-axis range. OO
filenames include centrality to distinguish plots from different selections.

| Macro/function | Input and purpose | Configurable arguments after the input |
|---|---|---|
| `plotEEChistogramsJetPtComparison` | Projected OO MC; compare the first four jet bins at one track threshold using the general drawer. Current cards give 40–80, 80–120, 120–160, and 160–200 GeV. | Track threshold, output directory, centrality minimum/maximum. |
| `plotEEChistogramsTrackPtComparison` | Projected OO MC; compare the first four track thresholds in the first jet bin. Current cards give 0.5, 1, 1.5, and 2 GeV in 40–80 GeV jets. | Centrality minimum/maximum; output goes to `figures/`. |
| `plotOOMCJetPtComparison` | Wrapper for `plotEEChistogramsJetPtComparison` with the same drawing style. | Output directory, track threshold, centrality minimum/maximum. |
| `compareOOMCPerpendicularBackground` | Projected OO MC with subevents and perpendicular-cone pairings; compare the estimated background with the same-jet mixed plus background-background components at 40–80 GeV and track pT above 2 GeV. | Centrality minimum/maximum; output goes to `figures/`. |
| `compareOOMCPerpendicularSignal` | Projected OO MC with subevents and perpendicular-cone pairings; compare total and background-subtracted EECs with the same-jet Pythia-Pythia component at track pT above 2 GeV. | Jet pT minimum/maximum, centrality minimum/maximum; output goes to `figures/`. |
| `compareFourOOMCModes` | Directory of projected RecoReco, RecoGen, GenReco, and GenGen files; compare total EEC shapes for track thresholds 0.5, 1, 2, and 3 GeV at 0–30% and 40–80 GeV. | Output directory. |
| `plotOOMCGenRecoTrackSelections` | Projected GenReco OO MC; draw the total EEC for the same four track thresholds at 0–30% and 40–80 GeV, without assuming reconstructed tracks have truth-origin labels. | Output directory. |
| `plotOOMCSubeventDecomposition` | Raw seven-dimensional RecoGen/GenGen OO MC EEC sparse histogram; project exponent 1 and plot Pythia-Pythia, Pythia-Hijing, and Hijing-Hijing components at 0–30% and 40–80 GeV. | Output directory, whether to use the verified stored generator range instead of resolving the four low-pT thresholds. |
| `plotRawOODataComparisons` | Projected OO data; compare four jet bins and four track thresholds for every stored centrality. | Output directory. |

Example calls for all nine macros are:

```bash
mkdir -p figures
root -l -b -q 'plotting/plotEEChistogramsJetPtComparison.C("projected/OO_MC_processed.root",2.,"figures",0.,30.)'
root -l -b -q 'plotting/plotEEChistogramsTrackPtComparison.C("projected/OO_MC_processed.root",0.,30.)'
root -l -b -q 'plotting/plotOOMCJetPtComparison.C("projected/OO_MC_processed.root","figures",2.,0.,30.)'
root -l -b -q 'plotting/compareOOMCPerpendicularBackground.C("projected/OO_MC_processed.root",0.,30.)'
root -l -b -q 'plotting/compareOOMCPerpendicularSignal.C("projected/OO_MC_processed.root",40.,80.,0.,30.)'
root -l -b -q 'plotting/compareFourOOMCModes.C("projected/fourModes","figures")'
root -l -b -q 'plotting/plotOOMCGenRecoTrackSelections.C("projected/OO_GenReco_processed.root","figures")'
root -l -b -q 'plotting/plotOOMCSubeventDecomposition.C("OO_RecoGen_raw.root","figures",false)'
root -l -b -q 'plotting/plotRawOODataComparisons.C("projected/OO_data_processed.root","figures")'
```

`compareFourOOMCModes` expects
`merged_OOMC_<mode>_20261022_processed.root` in its input directory, where
`<mode>` is `RecoReco`, `RecoGen`, `GenReco`, or `GenGen`. Its labels describe
the current generator-particle samples, which retain particles only from
5 GeV upwards; thresholds below that floor select the same stored population.
Use inputs matching that assumption. The raw subevent macro rejects inputs
that cannot resolve its requested low-pT cuts. For the existing generator skim,
pass `true` as its third argument to request a single diagnostic explicitly
labeled with the verified 5 GeV input floor:

```bash
root -l -b -q 'plotting/plotOOMCSubeventDecomposition.C("OO_RecoGen_raw.root","figures",true)'
```

The shape comparisons normalize each total EEC to unit area using the full
regular-bin width integral. Subevent components share the total normalization,
so their sum remains the total; the raw subevent macro checks both content and
variance closure before drawing. Projected mixed-subevent histograms already
combine both orderings, so use that component once. The perpendicular
background and signal comparisons instead retain the common raw EEC scale:
their estimate is the average of jet–perpendicular-cone pairs, minus pairs
between the two cones, plus the average of pairs within either cone. They do
not fit or independently renormalize the truth and estimate.

These plots diagnose shapes, available subevents, and background estimates.
They do not establish detector corrections or independent physics validation.
Reco-jet and gen-jet selections also contain different jet populations. Missing
truth-origin labels or particles below an input skim threshold cannot be
recovered by plotting or by changing the nominal track threshold.

### Framework support and compatibility

OO MC subevents are loaded from raw histograms, written during projection,
and available when reading projected files, for both individual jet bins and
the integrated jet selection. The drawer skips track processing when no track
plots are requested, allowing EEC-only inputs that do not contain the event
counter used by track QA. Track QA preserves variable EEC pT boundaries;
when grouping uniform pT bins it rescales the merged density so the width
integral remains unchanged. The drawer also initializes its maximum-particle
selection flags using their own array bound.

Shared sparse-projection helpers and track/particle matching readers use
`THnBase` so both float and double sparse inputs can be projected without
assuming `THnSparseD` storage. Matching projections also reserve space for the
optional energy-weight-axis restriction and apply it only when the histogram
contains that axis. The matching-closure getter keeps an already projected
histogram and reads from the file only when the cache is empty; this prevents
a valid raw projection being replaced by a missing processed-file object.
Multiplicity and particle-density arrays, initialization, and copying all use
the EEC track-bin capacity, preventing out-of-bounds access when an EEC card
has more thresholds than the ordinary track QA bins. Boundary-based EEC
lookups retain OO centrality selections and reject unavailable bin boundaries.

This plotting update preserves the explicit Run-3 response configuration,
native jet-bin selection, and collision-system/centrality compatibility checks
described below. Removing them would select legacy pp response files for OO
or ppRef, request jet-bin boundaries absent from the current cards, and drop
OO centrality selections. The plotting build retains `LinkDef.h` for the ROOT
dictionary and direct `TMath.h` includes in the model histogram managers so a
clean build has the required declarations. The existing
`plotting/validateOOMonteCarloWeights.py` and
`plotting/validatePpMonteCarloWeights.py` scripts remain available for marginal
weight checks; their same-sample one-dimensional comparisons do not establish
event-level closure.

### Validation of this plotting update

The reviewed commit was checked in an isolated source copy with ROOT 6.40.04:

- The analysis executable and a full rebuild of `DrawingClasses.so` compiled.
- ROOT loaded all nine new macros and the updated `plotEEChistograms.C`.
- `validation/checkOOCentralityLookup.cxx` passed, including invalid-boundary
  rejection and legacy/current histogram naming.
- `validation/run3Unfolding/testUnfoldConfiguration.cxx` passed, preserving
  explicit OO/ppRef response selection and legacy defaults.
- Synthetic QA checks passed with `MatchEECBinning` absent, disabled, and
  enabled, including axis edges, overflow storage, and preserved closure cuts.
- A synthetic projection/write/reload check passed 420,864 numerical content
  and error comparisons across float/double sparse storage, exponents 1/2,
  OO centralities, subevent sums, matching responses, and matching closure.

These checks verify software behavior; detector-response and production-MC
closure studies remain part of the validation work described below.

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

