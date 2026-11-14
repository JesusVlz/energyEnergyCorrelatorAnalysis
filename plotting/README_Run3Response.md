# OO Run-3 response diagnostics from projected campaign files

These additive C++/ROOT adaptations leave all Run-2 macros, the analyzer,
cards and projection code unchanged. They continue after jet-selection QA
and tracking comparisons. They do not perform unfolding or derive corrections.

## Run from repository root

```bash
make -C plotting
PROJ="${PWD}/projected/Run3MC_R4_20261112_v1_complete/OOMC" \
OUT="${PWD}/figures/OO_pairResponse_$(date +%Y%m%d_%H%M%S)" \
  bash runOOResponsePlots.sh pairs
```

Other families are `jets` and `unfoldingInputs`. Multiple names can be passed;
no arguments selects all three. Use a fresh OUT directory. Inputs are opened
read-only, logs and code provenance are recorded, and the driver stops on an
error or missing completion marker. No reprojection is run.

| Family | Configuration suffix | Projection used | Entry point |
|---|---|---|---|
| pairs | trackPairResolution | n1/n2, including -m | trackParticleMatchingPlotterRun3.C |
| jets | jetResponse | _processed.root, including -c | plotJetPtResponseMatrixRun3.C |
| unfoldingInputs | truthGuard5 | n1/n2, including -r -a -u | plotUnfoldingInputsRun3.C |

Input pattern: `${PROJ}/OO_MC_<scope>_<suffix>/<same-name>_n<N>_processed.root`.
Jet-only files omit `_n<N>`. Integrated scope selects the exact stored interval
[-1,100]; split selects 0-30, 30-50 and 50-90. No PbPb shift is applied.
Direct macros accept `(input, outputDirectory, centralityLow, centralityHigh)`;
pairs and unfolding inputs accept an optional fifth argument `trackThreshold=2`.
RecoReco is required for pairs and diagnostic jets; RecoGen for truthGuard5.
Physical boundaries, rather than fixed bin indices, select the reporting bins.

## Pair response

Adapted from `trackParticleMatchingPlotter.C`: nine figures per activity/power,
covering angle response, momentum-product response and momentum reco/gen
ratio distribution for jets 40-80, 80-120 and 120-160 GeV. The threshold is on
the lower RECONSTRUCTED track pT of the pair. These are matched pairs, not an
efficiency map or a comparison of all reconstructed/generated populations.
No Gaussian fit or resolution correction is derived in this first diagnostic.

The 2D display retains the legacy fixed-X normalization: X is reconstructed,
Y is generated, and each reconstructed slice sums to one over regular gen
bins. It represents P(gen bin | reco bin, matched, in-range), NOT P(reco|gen).
Empty slices remain blank, with ROOT masks and NaN in CSV. Stored errors are
scaled at fixed normalizer; normalization covariance is not calculated.
Flow is recorded and excluded from the conditional normalization. Log axes
omit nonpositive-edge bins from display only; these remain in ROOT/CSV.

Angle response has no energy-power axis and repeats in n1/n2 projected files.
Momentum axes are `(pT1*pT2)^n`, in GeV^(2n); the one-dimensional variable is
`[(pT1*pT2)_reco/(pT1*pT2)_gen]^n`. This matched-pair ratio is not a division
of two independent EEC histograms. Its display uses regular-bin unit sum and
fixed-normalizer marginal errors. Signed probability inputs are rejected.

## Jet response matrices

Adapted from `plotJetPtResponseMatrix.C`. The legacy -c matrix begins at
50 GeV and cannot supply an exact 40-80 slice. Its actual full matrix and
available 80-120/120-160 slices in both directions are drawn. An
`UNSUPPORTED_NATIVE_INTERVAL` message reports 40-80; nothing is relabelled
or interpolated. This does not invalidate the other projected outputs.

The native `oneDimensionalJetPtUnfoldingResponse` in truthGuard5 retains the
card's migration support and supplies exact 40-80 slices when populated.
Both families contain matched entries, not standalone efficiencies. Slices
have regular-bin unit sum; raw contents, errors and flow are preserved.
This does not implement the Gaussian JES/JER extraction from
`constructJetPtClosures.C` or derive scale factors.

## Unfolding inputs, not unfolding execution

Eleven figures per activity/power when all slices are populated: native jet
matrix, six matched jet slices, EEC-weighted flattened response, measured and
truth EEC marginals, and stored covariance input. The response/covariance maps
use bin indices for readability. Raw coordinates and JCard are saved; the
flattened coordinate is not a jet pT in GeV. All guard bins are retained.
The native jet matrix is independent of n; the EEC response is not.
Stored within-jet covariance is not certified full event or cross-file ratio
covariance.

Two CSV audits record total minus matched marginal contents, including the
other axis's flow in the matched sum. They require JetPtWeight=0 and assign
NO error to the difference because the samples overlap. Inspect residuals,
weights and flow before calling them misses/fakes. No unfolding closure,
iteration choice, efficiency correction or regularization follows from these
plots. Centrality-split does not mean independent train/validation halves.

## Outputs and software tests

Each figure writes PDF, PNG, ROOT and CSV with detached raw histograms,
display transformations, masks, JCard and provenance. Existing plot files
are not overwritten. A PASS marker means execution, not physics closure.

```bash
make -C plotting
root -l -b -q \
 'validation/run3Plotting/testRun3ResponsePlotters.C("/tmp/oo_response_fixture_unique")'
```

ROOT 6.40.04 built the plotting library and ran 45 synthetic panels. Tests
cover both powers, integrated/split indexing, all reporting jet bins,
conditional direction and error scaling, empty slices, signed conditional
input rejection, invalid modes/thresholds/centralities, and the legacy/native
40-80 distinction. PNG layouts were inspected. The November campaign files
were not executed here. Synthetic outputs are not physics results.
