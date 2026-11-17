# Histogram projection and plotting

This directory contains the histogram manager, drawing library, ROOT macros,
and Python plotting tools. Run commands from the repository root in a ROOT
or CMSSW environment; Python tools additionally require PyROOT.

```bash
make -C plotting
```

`EECHistogramManager` projects and loads analysis histograms. Its unfolding
covariance loader accepts both `THnSparseF` and `THnSparseD` and reports a
missing or incompatible covariance object instead of using an unchecked cast.
The producer's optional precision setting is described in the
[source README](../src/README.md).

## Run-3 jet selection, calibration, and veto plots

`plotRun3JetQA.C` reads merged analyzer files produced with `Run3JetQA 1`.
It plots cumulative and independent selections, calibration stages for a
fixed selected jet population, selected WTA eta/phi, and detector-veto maps
and event counts. Selection ratios share the same event exposure; they do
not receive independent/binomial error bars. Jet ID disabled in pp MC is
explicitly labelled, and inactive ID stages are omitted.

```bash
mkdir -p figures/run3_jets
root -l -b -q \
  'plotting/plotRun3JetQA.C("merged_OO.root","figures/run3_jets/OO",-1,100)'
root -l -b -q \
  'plotting/plotRun3JetQA.C("merged_pp.root","figures/run3_jets/pp",-1,0)'
```

Activity ranges are half-open: `[-1,100)` selects integrated OO and `[-1,0)`
selects pp. Selected eta/phi plots require the six exact jet-pT intervals
10-20, 20-30, 30-40, 40-80, 80-120, and 120-160 GeV from the new producer.
Veto occupancy ranges must align with the stored activity boundaries
`{-1,0,30,50,90,100}`. Older files without the new histograms report missing
inputs and skip those panels. Plots use weighted yields per bin width;
veto occupancy and event counts are unweighted. PNG/PDF figures, a
`<prefix>_numerical.root` file, and veto-event CSV accompany the plots.

`plotOOJetVetoMap.py` displays the nonzero support of the configured OO veto
map and compares it with the pinned payload. Its provenance JSON must provide
`inputs.integrated.path` and `inputs.integrated.sha256`; the merged input must
contain the configured map and `VetoMapFile`/`VetoMapKey` metadata.

```bash
python3 plotting/plotOOJetVetoMap.py \
  --provenance /path/to/provenance.json \
  --output-prefix figures/run3_jets/OO_veto
```

## All stored jet-pT intervals

```cpp
plotRun3JetKinematicsAllBins(input, prefix, jetLevel,
                            activityLow=-1, activityHigh=100);
```

`jetLevel` is `"Reco"` for RecoReco inputs or `"Gen"` for GenGen inputs.
This macro reads native `inclusiveJet` histograms for the v6 ordinary-jet
profile: minimum pT 10 GeV, maximum 5020 GeV, WTA axis, |eta| < 1.6,
`MatchJets 0`, `JetUncertainty 0`, and `JetPtWeight 0`. It requires the
11 card intervals with edges `{5,10,20,30,40,80,120,160,200,300,500,5020}`
and the native 5 GeV pT bins through 500 GeV.

It draws pT, eta, and phi separately for each interval and a pT overview,
preserving weighted counts and Sumw2. The 5-10 GeV interval is outside the
selection. The 500-5020 GeV interval aggregates native overflow; it cannot
resolve sub-bins or isolate the exact upper endpoint. Numerical ROOT and
coverage CSV files accompany the figures. Use `activityLow=-1` and
`activityHigh=0` for pp, and create the prefix's parent directory first.

## Tracking and pair-correction payload views

```cpp
plotRun3TrackingTables(input, outputDirectory, system, workingPoint);
plotRun3PairEfficiency(input, prefix, system, selection,
                      activity=0, jet=-1,
                      sampleLabel="Candidate; full training sample");
```

Both tools accept `system="OO"` or `"pp"`. The tracking plotter displays
existing efficiency, fake, secondary, and multiple-reconstruction tables:
OO uses eta/pT/hiBin tables and pp uses eta/pT tables. It preserves payload
bin contents and errors without deriving or averaging corrections.

The pair plotter reads an existing residual n=0 pair payload and its support
mask. It requires exact 6-8 and 3-4 GeV particle-pT cells, matching collision
system and track selection metadata, and valid activity/jet indices.
`jet=-1` selects jet-independent pairs. Unsupported cells are omitted;
horizontal bars represent bin widths and finite-training uncertainty is not
estimated. This quantity uses pair yields, without EEC momentum weights.
Create the prefix's parent directory before running the pair plotter.

## Response plots and existing guides

`Run3ResponsePlotTools.h` supports OO and pp reference MC and labels pp
samples with disabled Jet ID. Native jet-response slices extend over all
11 stored intervals when `MinJetPtCut` is 10 GeV; flattened response and
covariance displays mark jet-pT block boundaries while retaining guard bins.
Pair-response plots retain their existing 40-160 GeV reporting intervals.

See [response plotting](README_Run3Response.md),
[tracking and EEC plotting](README_Run3Tracking.md), and
[MC weights](README_Run3MCWeights.md) for the existing workflows and input
contracts. Generated plots, ROOT products, and study artifacts are separate
from the source files documented here.
