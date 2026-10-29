# Run-3 analyzer Jet ID and event veto

The user confirmed the pinned HIN PbPb Jet ID for OO2025 and pinned HIN pp Jet ID
for ppRef2024. Source commit: `c64293a132af388af1827044779af36e57f884d2`.
The exact cut table and differences from the earlier generic tight table are in
[vetoMapSources.md](vetoMapSources.md). The same functions apply to data and MC.

`Run3JetID 1` selects the collision-system HIN ID for reconstructed analysis,
reflected-cone QA and response jets. It replaces the inherited max-track/raw-pT
quality cut in these paths. Generator jets have no reconstructed PF ID applied.
Omitting the optional switch or setting it to zero preserves existing cards.
This disabled compatibility path does not certify the mandatory Run-3 veto.

`ApplyJetVeto 1` requires `Run3JetID 1` and all of:

| Card key | Required meaning |
| --- | --- |
| `JetVetoMapFile` | Explicit local ROOT path; no default or era inference |
| `JetVetoMapHistogram` | Explicit key; user nominated `jetvetomap_all` |
| `JetVetoPtConvention` | Exactly `raw`, `stored`, or `recomputed`; applicability must establish the choice |
| `JetVetoInputMinPt` | Declared complete input bound in that same pT convention, between 0 and 15 GeV |

The input-bound setting is an analyst declaration, not an automated proof that
upstream filtering retained every required jet. A known 20 GeV producer threshold
cannot satisfy it. The owner-enabled R=0.4 local and CRAB cards select stored
nominal PAT/HiForest `jtpt`, `JetVetoInputMinPt 15` and `jetvetomap_all`.
Both map sidecars record `ENABLED_FOR_PRODUCTION_TESTS`; full producer coverage,
map run coverage and exact campaign applicability remain unestablished.
The tarball script includes both maps and their metadata.

An event is rejected when ANY reconstructed jet has selected pT strictly above
15 GeV, passes its collision-system ID, has CEF+NEF strictly below 0.9, and lies in
a nonzero map bin. Negative map values veto too. The scan precedes EEC jet pT,
WTA eta, matching and pair selection and includes jets outside EEC acceptance.
The helper source's per-jet mask is deliberately aggregated into the requested
whole-event veto. This minimum selection does not add the legacy trackMax cut.

Map and ID use standard detector `jteta/jtphi`, while the EEC retains its selected
WTA axis. Recomputed veto pT, if explicitly chosen, preserves the existing selected
axis JEC lookup and nominal correction chain. Source evidence must approve that
combination. The loader requires explicit TH2 axes and finite bins; no underflow,
overflow, interpolation or edge acceptance fallback exists. Phi is periodic,
including Float_t rounding of plus/minus pi.

Ordinary readers bind native PF fractions/counts only when requested, preserving
original reco indices. Required types, dimensions and entry counts are checked;
nref/ngen bounds are read before payloads, including TChain transitions. Missing
inputs fail visibly. Charged counts are CHM+CEM+MUM, neutral counts NHM+NEM and
constituents their five-category sum, exactly as the selected HIN pp helper.
The generic-table reference API separately records HF-count incompleteness.

GenReco/GenGen signal events use an independently opened reconstructed companion
for the event veto. Mixing with generator tracks has its own synchronized reco
companion. The common event cut filters both signal and mixing with the same map
and predicate before pool matching. Mega skims cannot enable this veto because
they cannot establish complete eligible reconstructed jets.

Enabled outputs write unweighted signal/mixing checked, kept and vetoed counts
and eligible-jet detector eta/phi maps before and after whole-event rejection in
`eventQA/`. After maps exclude every jet in a rejected event. Disabled outputs
have no measured-loss histogram and explicitly record `NOT_APPLIED`. The card
persists the exact file, key and pT-convention strings under `JCard/*String`.

Reproduce software checks from the repository root:

```bash
python3 validation/run3Jets/run_selection_tests.py --output /tmp/NEW_EMPTY_TEST_DIRECTORY
```

The suite builds the analyzer and tests HIN cut boundaries, map validation,
reader schema/array failures, WTA/detector separation, all four sample types,
whole-event rejection below/outside EEC acceptance, signal/mixing consistency,
generator companions, configuration failures, saved QA and existing synthetic
response behavior. Private-access compilation is limited to regression fixtures.
These are software checks; they do not establish full-campaign detector closure,
map IOV, event-filter/trigger equivalence or a production selection efficiency.

## Executed Goal 3 checks

[Software checks](evidence/2026-10-28/software_checks.json),
[inclusive regression](evidence/2026-10-28/inclusive_regression.json) and
[native replay](evidence/2026-10-28/native-replay/native_replay_summary.json) are
retained. Replay uses the current C++ core and both explicit raw/stored pT
conventions. All retained jets were already above 15 GeV, so equal masks do not
resolve the convention or prove completeness. Rejected retained entries were
2/1000 OO data, 1/100 OO local MC, 10/100 pp data, 14/100 supplementary pp MC,
and 22/100 current pp MC. These are not analysis-selected event-loss estimates.

Two native charged fractions were one Float_t ULP above one. Validation accepts
two Float_t epsilon around physical endpoints without clipping or changing any
stored value used in an ID cut. Clearly invalid and nonfinite inputs still fail.

The owner confirmed keeping 15 GeV retention in ppRef data and MC. This is recorded separately from the inspected reference PSet's
20 GeV setting. Exact pT convention and other upstream retention cuts remain to verify;
a 15 GeV cut suffices only in the same veto convention with full eta coverage.

## Bounded MC executable validation, 29 October 2026

The current R=0.4 RecoReco executable passed 1,000-event checks on each local
OO and ppRef MC input. The retained results in
[`results/2026-10-29-mc-veto/`](results/2026-10-29-mc-veto/) record 962 selected
OO events, 386 jet entries and 15,448 EEC entries; ppRef has 843 selected events,
1,536 jet entries and 107,520 EEC entries. Signal veto counters are respectively
995 checked / 33 vetoed and 993 checked / 150 vetoed, after the preceding event
cuts. These counts are not campaign selection efficiencies or physics closure.

Use the [README commands](../../README.md) or
`python3 validation/run3Jets/run_mc_smoke.py --events 1000 --output NEW_DIRECTORY`
to reproduce with fresh outputs. The runner preserves the enabled map settings
from the base MC cards while bounding the event count and disabling backgrounds.
Historical manifests retain their original paths; the runner generates commands
for the current checkout and chosen output directory.
