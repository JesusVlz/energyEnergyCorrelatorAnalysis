# Run-3 R4 track and pair procedure

This directory applies AN2023_081_v13 Sec. 7.5.1 through the existing analyzer,
tracking tools, `TrackPairEfficiencyCorrector`, and `JDrawer` conventions.
Scope: OO selected MB (including stored hiBin 180–199), 0–30/30–50/50–90%,
ppRef, R4, reconstructed jet pT 40–80/80–120/120–160 GeV, charged tracks
pT>2 GeV and |eta|<2.4, momentum exponents 1 and 2. The OO activity
interpretation remains conditional on its established calibration.

## Native selections and one bounded comparison

`export_forest.py` runs native RecoReco and RecoGen on the same input entries.
Both use the analyzer's reconstructed-jet selection, WTA geometry, applicable
R4 JEC, HIN Jet ID, whole-event veto, event cuts and nominal track corrections.
The wrapper requires exact agreement of selected event identities, jets and
weights before merging. `--parity` compares the native EEC ROOT sums to exported
pair sums for both exponents. Source event/track IDs and full selected-event
track populations are retained. `export_forest_legacy.py` is preserved only to
reproduce the historical discrepancy; `compare_exporters.py` compares it to the
new export on identical bounded entries.

The established Run-3 analyzer uses the fitted vz weights, fitted OO centrality
weights within 0–90%, and its existing no-centrality-reweighting prescription for
90–100%. Its generator normalization factor is one; the stored forest weight is
retained as a separate field. This change does not establish a new cross-section
weight calibration. The former exporter used the stored generator weight and
omitted those fitted event weights; it also used stored jet pT and missed event
selection/Jet ID/veto differences.

Single-track factors now use physical detector eta and source-event hiBin
before any cone transformation. The installed OO tool retains its 6.1 GeV pT
lookup cap. For mixing, donor tracks keep their own factors/activity. Signal
self-pairs use the jet-dependent residual table; perpendicular/mixed-cone
self-pairs use the jet-independent table and source activity. Pairs formed
between distinct detector events or separated translated cones do not acquire
a close-pair correction from the artificial translated separation. This work
keeps the existing Run-2 correction path unchanged.

## AN residual maps and independent closure

`pair_counts.py` fills generator and single-track-corrected reconstructed pair
counts with event weights and **no EEC momentum weights**. The differential map
uses higher/lower track pT, source-event activity, reconstructed jet pT and
physical pair DeltaR. Both populations occupy the same reconstructed jets.
The AN jet-independent method uses full selected-event tracks, with its
higher-pT trigger inside |eta|<1.6 and associated track inside |eta|<2.4;
it is not an average of near-jet maps.

Within each differential cell, the control integrals use the strict, fixed
interval 0.1<DeltaR<0.4. If alpha = sum(reco_single)/sum(gen) there,
the multiplicative residual is C = alpha*gen/reco_single. ROOT tables use the
legacy normalized reco/gen histogram convention, so the existing consumer
applies its inverse. The minimum raw gen/reco counts is 10 in both cell and
control. Counts/control integrals and support masks are saved. This support
criterion is not a substitute for a finite-training-statistics uncertainty.
There is no smoothing, clamping, extrapolation or sparse-bin unity replacement.
Unsupported differential cells may use only a supported AN jet-independent
cell; otherwise the corrected value is unavailable. Invalid auxiliary full-event
single-track lookups are explicitly flagged and mask any affected independent
cell/control; an invalid lookup inside a selected jet stops the export.

The namespace/run/lumi/64-bit-event FNV1a64 modulo-2 split is unchanged:
fold 0 trains both tables, fold 1 validates. Physical duplicate events are
rejected across input files. Momentum factors `(pT1*pT2)^n`, each single-track
factor, residual factor and event weight appear once; detector corrections
are not raised to n. Results include raw, single, pair, jet-independent and
truth spectra, event vectors, table support and actual fallback use.

`plot_closure.py` writes compatible per-system correction ROOT files, closure
ROOT/PDF plots and representative AN angular/momentum-weight matching checks.
The latter reuse the established charge/dR/pT greedy matching on common selected
reconstructed jets. They do not establish negligible resolution effects or
replace the AN's subsequent smearing study. Closure statistical bars use common
physical-event contributions, conditional on the fitted table; finite-training
map uncertainty is not included. Unsupported whole-bin corrected yields remain
NaN; open points explicitly show only the supported partial contribution.

## Commands

From the repository root (ROOT/PyROOT required):

```bash
make -j4
python3 -m unittest discover -s validation/run3Pairs -p 'test_*.py' -v
python3 validation/run3Pairs/export_forest.py INPUT.root bounded/OO.jsonl \
  --system OO --sample-namespace CAMPAIGN_NAMESPACE --max-events 10000 --parity
python3 validation/run3Pairs/compare_exporters.py \
  bounded/OO_legacy.jsonl bounded/OO.jsonl bounded/OO_reconciliation.json
python3 validation/run3Pairs/run_sample.py RESULTS/expanded \
  --oo-files 34 --pp-files 10 --events-per-file 30000 --workers 4
python3 validation/run3Pairs/pair_counts.py RESULTS/expanded/OO_inputs.txt \
  RESULTS/expanded/OO_pairs.json --input-list --min-pairs 10
python3 validation/run3Pairs/pair_counts.py RESULTS/expanded/pp_inputs.txt \
  RESULTS/expanded/pp_pairs.json --input-list --min-pairs 10
python3 validation/run3Pairs/plot_closure.py RESULTS/expanded/OO_pairs.json \
  RESULTS/final/OO --events RESULTS/bounded/OO.jsonl
python3 validation/run3Pairs/plot_closure.py RESULTS/expanded/pp_pairs.json \
  RESULTS/final/pp --events RESULTS/bounded/pp.jsonl
python3 validation/run3Pairs/summarize_results.py \
  --oo-pairs RESULTS/expanded/OO_pairs.json --pp-pairs RESULTS/expanded/pp_pairs.json \
  --plots RESULTS/final --output RESULTS/final
```

The fixed local sample runner records commands, full input SHA256 hashes,
source/executable hashes and per-file receipts. `--resume` reuses an export only
when its command, source hashes, input bytes, output bytes and metadata still agree. `--input-list` streams JSONL files without building
a combined event file. The per-system manifests are subsets of `expanded/inputs.txt`.
Process systems separately to bound memory. For the compact resolution check, pass `--events` with
the representative bounded OO/pp JSONL files; the same physical fold is used.

Opt-in application in the normal analyzer uses `Run3TrackPairCorrectionFile`
and `DisableTrackPairEfficiencyCorrection 0`. The file must contain the
`Run3PairCorrectionVersion` marker and standard JCard/histogram names. The
Run-3 consumer fails on unsupported cells instead of silently accepting unity.
Default analysis cards stay disabled while delivered tables have unsupported
regions or observed nonclosure; generator-level correlations remain disabled.

Actual commands, first 40–80 plots, expanded tables/grid, source snapshots and
observed limitations are under `results/2026-11-02-an751-v1/`. See its
`RESULTS.txt` and the plan's Track pair efficiency Result/next entry.

The completed run contains 979,631 selected OO and 255,902 selected ppRef events.
All 30 channels are plotted. All six ppRef channels have complete holdout support,
with integrated pair/truth ratios 0.9932–1.0484 over [.008,.4); the smallest-angle
40–80 GeV bin still shows nonclosure. All 24 OO channels have incomplete support.
`final/channel_summary.csv` distinguishes full ratios from partial contributions;
`unsupported_cells.csv` and `control_support.csv` give exact cells and counts.
`audit_event_identity.py` verified all 1,320,000 raw input entries against the native
ledgers with no duplicate namespaced identities or mismatches.

## Preparing all local MC

`all_mc_campaign.py prepare OUTPUT` discovers actual filenames, reads full-file
ROOT entry counts with PyROOT and freezes sources, cards, calibrations and inputs.
It does not run MC. `run OUTPUT --workers 4 --resume --system all` executes the
frozen exports; `status OUTPUT` exits zero only when the full campaign is complete.
Successful receipts and native artifacts are verified before reuse; failures are
preserved. A new destination or changed source requires a fresh preparation.

The 3 November plan and staged commands are in the Track pair efficiency sheet
and `results/2026-11-03-an751-all-mc-v1/campaign_commands.sh`. Execution remains
pending: the measured home quota cannot hold the projected outputs, and the full
pp JSON/vector pipeline needs substantially more memory than the interactive node.
See `resource_readiness.json`; input streaming does not bound retained result memory.
