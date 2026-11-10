# Run-3 closure plots derived from the Run-2 macros

The Run-2 files `compareDifferentEECpairings.C` and
`validateBackgroundSubtraction.C` remain byte-for-byte unchanged. The Run-3
copies use the same ROOT/JDrawer structure, with explicit OO selections and
5.36 TeV metadata. No Python is required for generation, plotting or tests.

## Entry points

- `compareDifferentEECpairingsRun3.C`: raw-scale background, signal and
  method comparisons. Signal closure can be drawn with legacy independent
  histogram errors or a clearly labelled shared-signal approximation.
- `validateBackgroundSubtractionRun3.C`: the dedicated signal-validation
  copy. Its default retains the legacy independent unit-area normalization.
- `Run3ClosureTools.h`: input checks, detached histograms, numerical exports,
  provenance, and the explicitly limited shared-signal error calculation.

Only four superseded bespoke closure macros are removed:
`plotOOMCBackgroundPresentation.C`, `compareOOMCBackgroundMethods.C`,
`compareOOMCPerpendicularBackground.C`, `compareOOMCPerpendicularSignal.C`.
Other OO QA, subevent and response studies are retained. The new comparison
macro does not reproduce the removed six-family presentation atlas: it is a
focused closure/method-comparison replacement.

## Run the two requested error variants

From the repository root in a ROOT environment:

```bash
make -C plotting
root -l -b -q 'plotting/compareDifferentEECpairingsRun3.C("merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root","signal",false)'
root -l -b -q 'plotting/compareDifferentEECpairingsRun3.C("merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root","signal",true)'
```

Defaults: OO 0-30%, 80-120 GeV jets, track threshold 2 GeV, R=0.4,
common raw scale. The energy exponent is taken from the projected input.
The two variants have exactly the same central values, upper-panel errors,
axes and subtraction definition. Only the non-reference ratio errors change.

The full comparison signature is:

```cpp
compareDifferentEECpairingsRun3(input, plotKind, correlatedRatio,
  minJetPt, maxJetPt, minCentrality, maxCentrality, trackThreshold,
  outputDirectory, averageSameCone, normalizeShapes,
  ratioMinimum, ratioMaximum, extendRatioRange);
```

`plotKind` accepts `signal`, `background`, `signalWithTotal`, `methods`.
The signal comparison includes true signal, total minus mixed background,
and total minus perpendicular background. The background comparison uses
true background, mixed estimate and perpendicular estimate. The method
comparison uses perpendicular as its reference and mixed as its numerator.

Examples:

```bash
# Background closure with an explicit wider ratio range to show all deviations.
root -l -b -q 'plotting/compareDifferentEECpairingsRun3.C("merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root","background",false,40.,80.,0.,30.,2.,"figures/run3LegacyClosure",false,false,0.8,1.2,true)'

# A different native jet-pT interval.
root -l -b -q 'plotting/compareDifferentEECpairingsRun3.C("merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root","signal",true,120.,160.)'

# Dedicated shape-validation copy: independently unit-area-normalized by default.
root -l -b -q 'plotting/validateBackgroundSubtractionRun3.C("merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root")'
```

The dedicated validation signature is:

```cpp
validateBackgroundSubtractionRun3(input, backgroundMethod, correlatedRatio,
  normalizeShapes, minJetPt, maxJetPt, minCentrality, maxCentrality,
  trackThreshold, outputDirectory, averageSameCone);
```

`backgroundMethod=0` is mixed, `1` perpendicular, `-1` reads a previously
processed signal. Shared-signal errors require a reconstructed background
and `normalizeShapes=false`; they are rejected for the normalized-shape and
already-processed-signal modes because the needed covariance/provenance is
not supplied. Method `-1` is only a legacy comparison convenience: check the
input processing convention independently; the automatic background-definition
label is not a certification of its provenance.

## Preserve legacy arithmetic explicitly

The default `averageSameCone=false` reproduces what the current legacy macro
actually computes:

```
B = 0.5*(JC1 + JC2) - C1C2 + C1C1
```

The legacy source computes an averaged same-cone helper but adds C1C1 alone.
That behavior is not silently corrected in this first visual comparison.
Setting `averageSameCone=true` explicitly chooses:

```
B = 0.5*(JC1 + JC2 + C1C1 + C2C2) - C1C2
```

This changes finite-sample central values and marginal errors. It is separate
from choosing a correlated ratio error: never compare the two error options
with different `averageSameCone` values. Output names record `legacyC1` or
`averagedCones`. This is also why the default plots can differ from earlier
OO presentation plots that used the averaged definition.

No independent normalization is applied in the raw closure modes.
`normalizeShapes=true` is a separate, explicitly labelled shape comparison;
normalization and display ranges are not interchangeable.

## Ratio uncertainty model

The baseline is the existing ROOT `TH1::Divide` independent-histogram
propagation. The second signal variant uses:

```
T = S + Btrue
X = T - Best
R = X/S
Var(R) = [Var(Btrue) + Var(Best) + (R-1)^2 Var(S)] / S^2
```

The identical Pythia-Pythia signal is shared. Other cross-covariances are
set to zero in this approximation. It uses stored pair-fill Sumw2 marginal
variances, not a measured event-level covariance. Shared tracks, event weights,
background-component correlations, mixed-donor reuse and systematic errors
are not accounted for. The figure and ROOT metadata explicitly say
`shared-signal approximation (not full covariance)`.

`T=S0+S1+S3` and the marginal variance identity are checked before this
approximation is applied. Projected S1 already combines both mixed-origin
orderings; S2 is not added again. Unsupported background ratios are rejected
rather than assigned an invented correlated error. A full covariance-aware
closure remains a separate analysis task.

## Selections and output

Selections are resolved from JCard physical boundaries, not hard-coded PbPb
indices. No PbPb MC +4% shift is applied. For the supplied split input,
0-30% resolves to C1 because the first interval is -1 to 0. Exact unavailable
track thresholds and missing jet intervals fail explicitly. Truth-origin
closure requires OO RecoGen/GenGen with the required subevents. R=0.4 and a
single projected exponent (1 or 2) are required by these bounded copies.

Each plot writes PDF, PNG, ROOT canvas/histograms/ratios/provenance, and CSV.
Nonpositive denominators in the displayed range fail the requested figure;
they are not presented as zero measurements. The old fixed ratio zoom is
retained by default. `extendRatioRange=true` is an explicit display option.

## Reproduction and checks

`validation/run3Plotting/makeRun3ClosureCopies.cxx` is the checked C++
copy-and-patch recipe. Rebuild the copies without editing the originals:

```bash
g++ -std=c++17 -O2 validation/run3Plotting/makeRun3ClosureCopies.cxx -o /tmp/makeRun3ClosureCopies
/tmp/makeRun3ClosureCopies .
```

The recipe expects the source at commit
`792d76c0cd043ab1871b3217ac988f6455f2d9f0`; it fails on unexpected patch matches.
Preserved original blob hashes:

```
compareDifferentEECpairings.C  3ad7f8dfdb7eba8c936e8f9d9c78038ce262fb17
validateBackgroundSubtraction.C  c571a0ee20a72d9cc11cf7956d04b9a380476ffb
```

After producing both signal variants, run the C++ numerical test:

```bash
root -l -b -q 'validation/run3Plotting/checkRun3Closure.C("merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root","figures/run3LegacyClosure",80.,120.)'
```

Local validation used ROOT 6.40.04 and the supplied OO GenGen ROOT file.
The plotting library compiled, both error variants ran for 40-80, 80-120
and 120-160 GeV, and the numerical tests checked unchanged central values,
unchanged upper errors, the legacy background arithmetic, both ratio-error
formulas, and rejection of unavailable physical boundaries. The dedicated
shape-validation macro also ran. Legacy macros retain their existing Cling
variable-length-array warnings. These are software/plotting checks, not a
physics covariance validation. Input data and physics plots are not committed
to this public repository.
