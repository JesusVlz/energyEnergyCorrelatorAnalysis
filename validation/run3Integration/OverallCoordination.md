# OO/ppRef overall coordination: main integration

User-authorized software consolidation, 11 November 2026. This records the
repository work that is available, not inaccessible chat content or an assertion
that all analysis milestones are complete.

## Integrated work

The existing `oo-integration` lineage contains the tracked tracking/pair,
R4 unfolding, legacy-style Run-3 closure plotting, slim-input and exact PF-count
updates. Its base also contains the Run-3 R4 jet calibration/ID/veto/QA and
Nebraska storage support. Preserve their merge ancestry when promoting to main.
The R4 card generator is now accompanied by its generated catalog: 137 cards and
templates (75 configured, 62 requiring real inputs), plus two unfolding manifests.
See `cards/R4/README.md` and `catalog.json`; these are configurations, not jobs.

Additional changes in this integration:

* OO and ppRef at 5.36 TeV (data and MC; DataType 9/10/11/12) do not enable,
  bind or consume `trkNormChi2`, `trkNLayers`, `trkNHits`, `pfEcal`, `pfHcal`.
  Their AOD chi2/ndof/hit/layer equivalents are also skipped. This applies to
  standard and mixed-event HighForestReader paths and UnfoldingForestReader.
  Accidental getter calls fail explicitly, even if a full forest has those fields.
* Active Run-3 cards omit `CalorimeterSignalLimitPt`, `HighPtEtFraction`,
  `Chi2QualityCut`, `MinimumTrackHits`. The analyzer does not request these keys
  for Run-3. Generators cannot reintroduce them through an older base card.
  Run-2 configurations, archived snapshots and historical result/evidence cards
  retain their original settings. This is not a change to the track prescription.
* `plotting/deriveMonteCarloWeights.C` calls the new Run-3 implementation.
  The original Run-2 implementation remains callable as
  `deriveMonteCarloWeightsRun2.C`. Stored-weight plots and candidate 1D refits
  use JDrawer; no fits or correction parameters are deployed automatically.

## Preserved analysis contracts

R4 WTA, integrated and split OO, ppRef, reporting at least 40-80/80-120/120-160
GeV, both momentum powers, and migration-support bins remain distinct.
Integrated OO uses the ordinary centrality axis/projection, not InclusiveMode.
R8 remains deferred outside the active catalog; no AK4-to-AK8 substitution.
Jet ID is additional to trackMax/rawPt. pp L1FastJet stays off for existing
Vanderbilt forests without the required native inputs. Signal and donor lists
are chosen independently: Nebraska mixing does not relocate signal files.

Single-track corrections use physical detector coordinates and source-event
activity with the matched Nominal/Loose/Tight payload. Pair maps are derived
from unweighted momentum-pair counts in common reconstructed jets; training
and holdout remain separate. Keep unqualified pair corrections disabled.

The sequence remains track/pair-corrected EEC -> jet-pT unfolding -> adjusted
background subtraction -> normalization -> OO/pp. Compare mixed and perpendicular
backgrounds with explicit exposures. Existing pp no-nominal-subtraction policy
is retained; pp MC background checks are diagnostic, not assumed zero background.
Response preparation and MC closure do not complete measured-data systematics.

## Weight bookkeeping

OO already applies vz and centrality fits; ppRef already applies its own vz fit.
Run-3 generator/cross-section factor remains unity in the event loop.
`vertexZ`/`centrality` carry that baseline factor; `vertexZweighted` and
`centralityWeighted` carry the full event weight. Never multiply candidate
replacement fits into already fully reweighted marginals. Stored old files
retain the weights of the executable that produced them. A 1D same-sample
refit is not independent closure of the combined OO product.

## Evidence and remaining execution

`bash validation/run3Integration/checkMainUpdate.sh` builds the analyzer and
plotting library; checks all active/generated cards, synthetic slim/full
AOD/MiniAOD reader inputs including chained mixed events, table/selection/jet/PF
regressions, unfolding software tests, closure macro compilation, and synthetic
ROOT weight fits/figures. It never submits production. Inspect the workflow
result and artifacts for the exact tested SHA; missing RooUnfold is an explicit
runtime skip, not a passed covariance test.

Native OO/pp MC runs, independent correction/response/full-chain closure,
qualified pair payloads, data fits and measured-data systematics remain actual
execution tasks. Preserve hashes and ledgers and do not relabel old products.

After syncing the tested main revision, rebuild in the compatible CMSSW/ROOT
environment (`make clean && make -j2; make -C plotting -j2`). Existing binaries,
CRAB sandboxes and campaign tarballs are not changed by a Git update; rebuild
and repackage those before using the slim forests. Do not add unused branches
back merely to silence an older executable.
