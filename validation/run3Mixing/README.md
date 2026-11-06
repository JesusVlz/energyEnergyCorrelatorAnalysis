# Run-3 two-partner mixing preflight

MC background item-list handoff for ChatGPT and continuation: [MC_BACKGROUND_HANDOFF_2026-10-29.txt](MC_BACKGROUND_HANDOFF_2026-10-29.txt).

`preflight.py` measures pool feasibility using native jets and stored measured
`hiHF_pf` and vertex coordinates. It never fills or subtracts background EECs.
The HF/vz tolerances are explicit **study settings**, not an approved matching
prescription. It applies no HF shift, centrality replacement, matching weights,
pair correction, component normalization or clipping.

The bounded scan saves every selected event, including zero-jet events, to an
adjacent `.events.jsonl` ledger. It retains real hiBin (including unsupported
values), signed generator weights, namespace/run/lumi/event identity, source
entry, native radius and jet bin. Source/accepted 0–90, 90–100 and invalid hiBin
counts are recorded independently of the matching coordinate. Invalid hiBin in
pp is bookkeeping, not an OO percentile acceptance diagnosis.

Required MC filters follow the current Run-3 reader: OO primary vertex,
`OOphfCoincFilterPF2Th4`, pileup; pp primary vertex and `ppNoScrapingFilter`;
abs(vz)<=15. Missing required branches fail. No data trigger is inferred. Source
skim losses and MB representativeness cannot be established from selected MC.
Native pp R8 absence is explicit and never replaced by AK4 jets. Stored jtpt is
provisional while applicable calibration remains unresolved.

Each signal physical event gets two distinct other physical events in the same
sample namespace and FNV training/holdout fold. Vertex and HF mismatches must
satisfy both requested tolerances; sparse pools fail without relaxing them.
The stored HF range and distinct-value count are audited. A constant or empty
coordinate is explicitly BLOCKED_DEGENERATE_COORDINATE; feasible partner counts
then describe vertex-only feasibility and cannot establish activity matching.
Candidates are ranked by normalized L1 mismatch, then stable event identity.
Partners are shared between radii/powers of the same signal event. Reuse counts,
failed signals, jet-bin coverage and signed weights remain visible. These records
support later covariance work; they are not a covariance calculation. In
particular, partner reuse is not independent exposure.

Example (choose a new output path; existing outputs are protected):

```bash
python3 -m unittest discover -s validation/run3Mixing -p 'test_*.py' -v
python3 validation/run3Mixing/preflight.py \
  /eos/user/j/jvelazqu/Foresting/OOMC/unembedded_L1HLT_objects_MC_1.root \
  /tmp/OO_pool_study.json --system OO \
  --sample-namespace OO_local_unembedded_L1HLT_MC_1 \
  --events 10000 --vz-tolerance 1 --hf-tolerance 20
```

Use the same verified physical-sample namespace as the pair exporter. The
example namespace identifies the local audited file, not certified campaign
provenance. Study tolerance sensitivity by importing `match` and rerunning it on
the saved ledger. Do not optimize tolerances on an OO/pp physics ratio.

Promotion requires independent partner-source applicability, matching and jet-HF
bias studies, adequate channel coverage, measured data/MC compatibility, component
exposures/correlations, R8 overlap treatment and extracted-signal MC closure.
The 2 November method comparison and the inherited final correction order remain
unchanged. No production campaign or nominal estimator is authorized by this
preflight.

## 2026-10-29: MC-first OO R4 GenGen closure pilot

The existing pure-HIJING list is
`mixingFileList/OO2025_HIJING_MC_pool_pilot.txt` (33 accessible native files).
The source/schema/provenance audit is
`OO2025_HIJING_MC_pool_pilot.manifest.json` beside it. The `_A.txt`/`_B.txt`
files are file references only; the authoritative partition is **event-level**,
using nonempty stored-background fingerprints, or full physical identity for
empty backgrounds. Different paths do not establish independence.

`cards/OO_GenGen_R04_closure.input` enables an opt-in path in the existing
analyzer. It uses native R4 WTA GenGen jets, preserves generator momenta,
asserts unity reconstruction track/pair corrections, and retains the existing
Run-3 signal selection and reconstructed-jet whole-event veto. The Run-3 jet
and tracking payloads remain in the applicable reconstructed paths.

Matching uses stable charged HIJING particles with strict pT>0.7 GeV and
abs(eta)<2.4, independent of the EEC strict pT>2 GeV cut; abs(dz)<=0.5 cm and
abs(dNH)<=max(1,0.01*NHsignal). These are pilot tolerances. HF transformations
are OFF with exact raw-HF preservation. The capped signal prescan globally
excludes donor raw event identities and observed nonempty signal-background
fingerprints. A/B donors must be distinct; finite searches fail without
widening. One assignment per physical signal event supplies every jet/power.
Pure-HIJING particles are background by source, despite their stored sube=0.

The paired estimators and total/PP/PH/HH share the same matched signal jets
and event weights. No donor hard-process weight is multiplied into pairs.
The raw manager API bypasses background-shape stabilization; existing manager
loading and reconstructed results retain their behavior. PH merges both origin
orderings exactly once. No independent area normalization or negative-bin
clipping is used. Conditional fixed-pool errors retain joint signal-event
vectors and observed background clusters; donor-resampling and final covariance
remain subsequent work. The current reconstructed pair routine uses relative
eta and signal hiBin; this pilot does not alter that routine's physics.

Run from the configured ROOT environment (choose a new output directory):

```bash
make -j4
python3 -m unittest discover -s validation/run3Mixing -p 'test_*.py' -v
python3 validation/run3Mixing/run_closure_pilot.py \
  --signal /eos/user/j/jvelazqu/Foresting/OOMC/unembedded_L1HLT_objects_MC_1.root \
  --output-dir validation/run3Mixing/results/2026-10-30-gengen-pilot/bounded_100k \
  --signal-events 10000 --pool-events 100000
python3 validation/run3Mixing/plot_closure_pilot.py \
  --input validation/run3Mixing/results/2026-10-30-gengen-pilot/bounded_100k/pilot.root \
  --ledger-prefix validation/run3Mixing/results/2026-10-30-gengen-pilot/bounded_100k/pilot \
  --output-dir validation/run3Mixing/results/2026-10-30-gengen-pilot/bounded_100k/plots
```

The runner saves cards, exact commands, caps/seed, source hashes/snapshot,
software versions, and native run logs. Results are local artifacts ignored by
Git. The ledger checker separately verifies physical IDs, source coordinates,
matching, exposure and reuse against the audit. Projection/plotting covers all
three requested jet bins, both powers and integrated/three split selections,
with legacy 0.008-0.39 and full stored R4 angular views. Empty channels remain
explicit; ratio masks retain zero truth denominators.

The exact planning references are `/eos/user/j/jvelazqu/HIN/Plan/To do.md`
and `/eos/user/j/jvelazqu/HIN/Plan/Mixed and perpendicular cones.md`.
AN2023_081_v13 and Jussi's 2025-07-14 slides were not found locally and are not
claimed as read. The original generator-background lineage and historical
producer retention settings are not fully pinned by the ROOT metadata; observed
stored-particle exclusions cannot certify latent lineage independence. Actual
selected hiBin acceptance is reported separately from a physical 0-100% claim.

This is the 30 October presentation pilot; the extended comparison milestone
remains 2 November. Measured-activity matching, reconstructed closure and final
covariance follow separately. The final sequence remains tracking/pair-corrected
EEC -> jet-pT unfolding -> appropriately adjusted background subtraction ->
normalization -> OO/pp. The estimator is not frozen.

### Completed evidence and submitted OO MC campaign

The authoritative bounded run is
`results/2026-10-30-gengen-pilot/bounded_100k/`: native execution, independent
100k donor/source comparison, full partner/exposure ledger checks, and raw
manager projections all passed. It matched 9,486/9,596 selected signal events;
the requested jet bins have 193, 9, and 0 matched jets. The 193 PDF/PNG plot
pairs are in `bounded_100k/plots/`, including explicit T/PP/PH/HH decomposition.
HH truth is zero throughout and 80-120 GeV truth background is zero; the empty
120-160 GeV bin and invalid ratios stay masked. Software validation passed;
statistics and conditional fixed-pool errors do not establish final closure.

Final checks: 30 mixing tests, 18 Condor tests, Run-3 jet regressions, all nine
forced-build workflow stages, and independent native identity/weight/exposure
checks. See `checks/final_checks.json` and `checks/run3Workflow/summary.json`.
The initial first-file-count attempt was stopped and retained; the additive
pilot total-chain API and two-file regression fixed that counting issue.

On the owner's request, separate full-file OO-only GenGen campaigns were
prepared and submitted under `/afs/cern.ch/user/j/jvelazqu/condor/production/OOMC/`:
`GenGen_R04_Run3_NH_20261029_split` (cluster 9624972) and
`GenGen_R04_Run3_NH_20261029_integrated` (cluster 9624973), 333 jobs each.
Each reads the full forest (cap 30,000), using the same fixed 100k donor prefix
and `OOClosurePartitionPolicy 1`: reserve all donor identities/backgrounds,
then exclude overlapping signals before common exposure. This campaign-wide
rule is distinct from the test's default policy 0 (reserve its capped signals).
Both policies keep the same NH/vz matching, strict EEC pT>2 and both powers.
Read `condor/README.md` for the frozen payload, local pool staging, resources,
publication completion marker, and preparation/worker tests. Never combine
split and integrated productions as independent exposure. The test card and
its 10k limit remain separate from these full-file campaign cards.

Exact submission records, frozen provenance and queue snapshots are in each
campaign. `results/2026-10-30-gengen-pilot/campaign_summary.json` records both
clusters/settings and the local evidence. The background task sheet has been
updated with actual code/run/physics status and the 2 November follow-up gates.
