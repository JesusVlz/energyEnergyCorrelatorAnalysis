# OO Run-3 tracking and EEC plotting

`trackClosurePlotterRun3.C` adapts the generated/corrected/uncorrected track
comparison, per-event normalization and JDrawer split-panel layout from
`trackClosurePlotter.C`. The original Run-2 macro is unchanged.

## Build and run

Run from the repository root in your ROOT/CMSSW environment:

```bash
make -C plotting
PROJ="${PWD}/projected/Run3MC_R4_20261112_v1_complete/OOMC" \
  bash runOOTrackingPlots.sh Nominal
```

The driver accepts `Nominal`, `Loose`, and `Tight` and reads:

```text
${PROJ}/OO_<workingPoint>_RecoReco/OO_<workingPoint>_RecoReco_n<N>_processed.root
${PROJ}/OO_<workingPoint>_RecoGen/OO_<workingPoint>_RecoGen_n<N>_processed.root
```

Both powers `N=1,2` are required. It requests centralities 0-30, 30-50 and
50-90, and EEC jet-pT bins 40-80, 80-120 and 120-160 GeV. Physical boundaries
are resolved from JCard, not assumed indices. This driver targets the paired
tracking-preparation configurations, not every configuration in the campaign.

A fresh default output directory is created at
`figures/OO_MC_tracking_<timestamp>/`. Set `OUT` to a nonexistent directory to
choose another location. The driver stops on errors, captures individual logs,
records the code version, and requires a successful plotting completion marker.

```bash
PROJ="${PWD}/projected/Run3MC_R4_20261112_v1_complete/OOMC" \
OUT="${PWD}/figures/OO_tracking_LooseTight_$(date +%Y%m%d_%H%M%S)" \
  bash runOOTrackingPlots.sh Loose Tight
```

There are 27 panels per working point: three centralities times three
inclusive-track panels plus three EEC jet-pT panels for each of two powers.
Inclusive-track panels are drawn only from n=1 to avoid duplication.
Each panel writes PDF, PNG, ROOT and CSV. Existing ROOT plot outputs are not
overwritten. Inputs are opened read-only; no reanalysis or projection is run.

## Direct first call

```bash
make -C plotting
PROJ="${PWD}/projected/Run3MC_R4_20261112_v1_complete/OOMC"
RECO="${PROJ}/OO_Nominal_RecoReco/OO_Nominal_RecoReco_n1_processed.root"
GEN="${PROJ}/OO_Nominal_RecoGen/OO_Nominal_RecoGen_n1_processed.root"
root -l -b -q \
  "plotting/trackClosurePlotterRun3.C(\"${RECO}\",\"${GEN}\",\"figures/OO_Nominal_tracking_first\",0.,30.,2.,\"Nominal\",true,true)"
```

Signature:

```cpp
trackClosurePlotterRun3(recoFile, genFile, outputDirectory,
  centralityLow=0, centralityHigh=30, trackThreshold=2,
  workingPoint="Nominal", drawTracks=true, drawEEC=true,
  angularTrackHigh=3);
```

## Plot definitions and uncertainty model

* Inclusive single tracks: generated particles, corrected tracks and
  uncorrected tracks, normalized per selected weighted event using the legacy
  integral-with-bin-width convention for multiplicity. Nonempty multiplicity
  flow is rejected. Track pT is displayed from 2 to 100 GeV. Eta and phi use
  the exact stored 2-3 GeV slice by default, not an integrated pT>2 selection.
  These are inclusive tracks in selected events, not tracks restricted to
  the EEC reporting jet bins.
* EEC: all-subevent, same-jet totals in RecoReco divided by RecoGen, with
  reconstructed jets in both. No independent area rescaling, background
  subtraction or unfolding is performed. Residual pair correction must be
  OFF. This is not extracted-signal/Pythia-Pythia closure and is not a method
  for deriving a residual pair-efficiency map. The display range is
  0.008 < Delta r < 0.39.
* EECs are not rebinned. Track pT and angles use grouping factors 10 and 2,
  implemented through bin-integrated sums to preserve densities. Output
  histograms are detached copies; inputs are not modified.
* Ratio errors propagate stored marginal histogram errors with zero
  covariance and fixed per-event normalizers. Figures explicitly label this
  independent-histogram approximation. No event-level or cross-file
  covariance is inferred, and no shared-signal background-closure formula
  is transplanted into these tracking ratios.
* Nonpositive denominators are omitted from graphs and recorded as NaN
  with valid=0 in CSV, not zero measurements. Ratio ranges expand to include
  displayed errors rather than hide deviations behind a fixed zoom.

## Input checks and limitations

Checks cover OO MC, RecoReco/RecoGen argument order, one projected exponent,
matching selected card entries, exact bins, R=0.4, all subevents, working-point
vertex cut, residual pair correction OFF, finite histograms, equal weighted
selected-event sums and compatible reconstructed-jet spectra.

The aggregate exposure/jet checks use relative tolerance 1e-4 for accumulation
order. They do not establish physical event-identity agreement, matching
calibration payloads or complete production coverage. Not every Run-3 setting
is persisted in the legacy JCard; inspect the corresponding merge receipts.
A `RUN3_TRACKING_PLOTTING_PASS` marker means software completion, not a
physics-closure verdict. Failures print `RUN3_TRACKING_PLOTTING_FAIL` and exit
ROOT with nonzero status.

## Synthetic software tests

Reviewed source interfaces: commit 6f7464a871d6c39e1b5e4211ed2fa01bf87c6a9e.
The macro and Bash driver are the versions supplied in OO_Run3_NextPlots.zip.
The test was relocated under validation/run3Plotting with its relative
include adjusted; no generated ROOT files or physics plots are committed.

The package's recorded ROOT 6.40.04 fixture run checked both energy powers,
all three reporting jet bins, the ratio/error formula, zero-denominator
omission and rejection of incompatible modes, powers and jet populations.
These are synthetic software checks, not results from the November campaign.

To rerun, use a fresh output path:

```bash
make -C plotting
root -l -b -q \
  'validation/run3Plotting/testTrackClosurePlotterRun3.C("/tmp/oo_tracking_fixture_new")'
```

The expected final marker is `SYNTHETIC_TRACKING_TESTS_PASS`.
