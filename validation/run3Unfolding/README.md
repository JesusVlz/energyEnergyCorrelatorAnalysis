# Run-3 R4 unfolding preparation

MC-only preflight, independent-fold scan and diagnostic plotting, based on
HIN-23-004 / AN2023_081_v13 Sections 7.5.2, 7.5.3 and 7.5.8. This is a preparation,
not a completed physics validation or frozen response.

## Integration

The entrypoint reuses EECCard, EECHistogramManager, EECUnfoldConfiguration and
RooUnfoldBayes. It does not implement a new unfolding estimator, fill responses,
resplit forests, submit jobs or change native selections/calibrations. Reuse the
existing batch/response_worker.py and batch/merge_response.py products.

Reviewed baseline: 792d76c0cd043ab1871b3217ac988f6455f2d9f0. Native files and Run-2
macros remain unchanged. **Use this new Run-3 entrypoint**, not the Run-2-default
iteration/plotting macros. The old GetMeasuredCov output path in
plotting/unfoldEEChistograms.C is NOT patched in-place. This adapter bypasses it
and retrieves the actual unfolded covariance through Eunfold or Ereco. The
installed API is recorded; an unsupported API fails, never falls back to input
covariance. ROOT integration still needs to run on the analysis host.

Scope: OO and ppRef, R=0.4, WTA, |eta_jet|<1.6, particles pT>2 GeV and |eta|<2.4,
both energy powers. Report 40-80, 80-120 and 120-160 GeV; retain native guard bins.
One manifest entry selects one existing activity interval and one projected power.
For split OO files use [0,30], [30,50], [50,90]. The batch selected-MB inclusive
projection uses [-1,100], including bookkeeping; it is not relabeled strict
0-100%. ppRef uses its native inclusive [-1,100]. No HF shift or R8 changes.

## Files

- run3_study.py: check or scan; explicit inputs, two fold directions, 1-20 iterations.
- study_root.py: native histogram/card interface and result-covariance adapter.
- study_math.py: support, conditioning, normalization and diagnostic statistics.
- draw_run3_study.py: response heatmaps, iteration curves, three-bin shape/yield closure.
- test_study*.py: numerical, metadata, API, ledger, renderer and ROOT regressions.
- examples/: input schemas, not actual production evidence.

## Prepare inputs

Use compatible **per-power projected train/validation** products from the existing
batch merger. Copy examples/study_manifest.json to a local study_inputs.json and
replace the paths. Add entries for available systems/activity/powers with unique
IDs. Paths may be absolute or relative to the input manifest.

`check` requires projected products but not split evidence. `scan` also requires
an index with the examples/split_evidence.json schema: hashes of the exact
projected products plus hashes and paths of their original worker identity
ledgers. Obtain hashes with sha256sum. This is a new study-side index, not a claim
that existing batch receipts already use that schema. Retain the original merge
and worker receipts as the evidence linking these products to those ledgers.
The checker verifies supplied hashes and global identity uniqueness across the
ledgers; it cannot independently establish hidden embedding-source reuse or
prove an incorrectly asserted product-to-ledger linkage.

Default truth_definition is **response_marginal**. For the separate ordinary
RecoGen-to-GenGen test set truth_definition to ordinary_GenGen and provide all
four train_reco, train_gen, validation_reco, validation_gen paths. Include those
four hashes in split_evidence.product_sha256 as well. All reference products
must use the same original event folds and compatible metadata. No response
marginal is substituted for a missing ordinary GenGen distribution.

Measured covariance is mandatory and must come from compatible RecoGen filling.
The native projection uses bin-integral second moments (not bin-width densities).
Its equivalence to the selected measured observable must be checked on actual
inputs. Ordinary EEC densities are converted to bin integrals before unfolding.
Do not pass previously normalized shapes as ordinary references.

Optional fake_scales with both A_to_B and B_to_A entries enable a fixed-exposure
refolding diagnostic. Establish the scales from the compatible weighted MC
exposures; do not tune them to obtain closure. Validate that convention against
the installed RooUnfold fake handling. The fake term here is unmatched jets, NOT
the track-pair background subtracted later. Without scales refolding is NOT_RUN,
not silently approximated with unit-column normalization.

## Run on the established ROOT/CMSSW analysis host

Requirements: current DrawingClasses.so built by the normal repository procedure,
RooUnfold with PyROOT dictionaries, Python >=3.8, NumPy and Matplotlib. This change
does not install software or rebuild the framework automatically.

```sh
python3 -m unittest discover -s validation/run3Unfolding -p 'test_study*.py' -v

# Require the real RooUnfold test on the ROOT host; missing ROOT must fail here.
EEC_REQUIRE_ROOT_TESTS=1 python3 -m unittest discover \
  -s validation/run3Unfolding -p 'test_study_root_runtime.py' -v

python3 validation/run3Unfolding/run3_study.py check \
  --manifest study_inputs.json --output validation/run3Unfolding/results/my-preflight
python3 validation/run3Unfolding/draw_run3_study.py \
  validation/run3Unfolding/results/my-preflight

python3 validation/run3Unfolding/run3_study.py scan \
  --manifest study_inputs.json --output validation/run3Unfolding/results/my-scan
python3 validation/run3Unfolding/draw_run3_study.py \
  validation/run3Unfolding/results/my-scan
```

Use fresh output directories. `--drawing-library` and `--roounfold-library` can
select existing libraries. Failure in one channel does not prevent other
channels being attempted; status is recorded and the overall exit is nonzero.
There is no collision-data unfolding or response-freezing command.

## Outputs and safeguards

Each checked channel writes manifest.json and plot_data.npz. A successful scan
also writes unfolding_scan.root, chi2_scan.csv, candidate_iterations.json. The
ROOT file retains input histograms, measured covariance and all 40 unfolded
results with actual truth-space covariance. Failed scans retain a clearly named
partial ROOT file, not a completed result. The top-level manifest records code,
input hashes, Git state, commands, errors and per-channel status.

Plotting a preflight draws only the two response matrices. Plotting a successful
scan adds split-1/split-2/average iteration curves and three-panel shape and raw
yield closure plots, saved as PNG/PDF. Raw yields remain available; shape
normalization cannot conceal yield failures. No placeholder physics plots are
created. Numerical and rendering tests use synthetic arrays, not MC evidence.

Native support and numerical rank are reported for K=C/T and conditional Kprime.
The adapter does not delete unidentifiable truth directions or clear flow bins.
Its current flow support is restricted to explicit regular guard bins; nonzero
ROOT flow blocks this entrypoint pending explicit modeling. That is an adapter
constraint, NOT a claim that ROOT flow bins are physically forbidden.

The default diagnostic DeltaR interval [0.008,0.39] follows the inherited
FindBin endpoint convention. The actual full-bin range is recorded; this is not
a frozen OO physics range. The scan uses the same valid-bin mask across both
folds and all iterations, sums the diagonal fixed-unity chi2 over all three jet
bins and records one diagnostic minimum per response. There is no inherited
four-iteration choice or separate optimum per displayed jet bin.

## Statistical limitations and CMS checklist

The measured per-jet covariance is propagated by RooUnfold; the new adapter gets
the actual unfolded covariance. Shape covariance uses the normalization Jacobian.
But truth errors currently use a diagonal Sumw2 proxy, and joint unfolded/truth
cross covariance is not available. Thus diagonal chi2 curves and ratio errors
are explicitly **diagnostic**, not calibrated p-values, coverage or final
iteration validation. Response-MC statistical uncertainty is not included in
these curves. The A/B tuning tests are not untouched final validation samples.

Still pending: actual ROOT/native-library integration and MC execution; input
retention and guard/binning stability; generator normalization and calibration
applicability; joint closure/refolding covariance and response-MC statistics;
unregularized benchmark where identifiable; prior/model stress, formal CMS
bottom-line test and coverage; alternative-generator assessment when available;
tracking/pair correction, adjusted post-unfolding background subtraction and
full-chain closure. The measurement order remains corrected EEC -> jet-pT
unfolding -> adjusted background subtraction -> normalization -> OO/pp.
