# Goal 3: Run 3 veto-map candidate evidence

Audit date: 2026-10-28. This task investigates the two maps explicitly nominated by the user. Source retrievals were kept under `/tmp/run3-goal3-veto-audit/`. The user subsequently uploaded byte-identical copies into `jetVetoMaps/OO2025/` and `jetVetoMaps/ppRef2024/`; their source hashes were verified and candidate-status metadata sidecars added. Neither map is enabled as nominal or included in a job package. The complete machine-readable record, including every key, full axis edges, file hashes, and unresolved fields, is [vetoMapCandidates.json](vetoMapCandidates.json).

## Sources and their limits

The prior calibration manifest pins [HIN-Jets at c64293a132af388af1827044779af36e57f884d2](https://github.com/thenicholasbarnett/HIN-Jets/tree/c64293a132af388af1827044779af36e57f884d2). The [pp helper](https://github.com/thenicholasbarnett/HIN-Jets/blob/c64293a132af388af1827044779af36e57f884d2/include/JetSelection_pp.h) and [PbPb helper](https://github.com/thenicholasbarnett/HIN-Jets/blob/c64293a132af388af1827044779af36e57f884d2/include/JetSelection_PbPb.h) load `jetvetomap_all` and evaluate eta on X, phi on Y. They combine a per-jet map decision with ID, rather than the requested whole-event veto. Their default points to `Winter25Prompt25_RunCDE.root`; that default supplies no OO2025/ppRef2024 approval and was not adopted. The pp ID thresholds also differ from the user-provided table; the PbPb helper applies only leptonic terms. The initial generic CHS table and these helper policies have different cuts. The user explicitly selected the linked collision-system helper policies; those are the analyzer implementation. The generic-table helpers remain isolated core regression references.

The latest user steering explicitly requests the PbPb-helper strategy for OO and the pp-helper strategy for pp, and states that the OO strategy is officially approved. That statement is recorded as owner-provided approval of the strategy. It does not supply a numerical map IOV, exact producer attribution, or a retrieved CMS approval record for a particular binary. A fresh lookup confirms current HIN-Jets `main` is still the pinned commit above. Neither helper accepts a pT argument, applies the separate CEF+NEF<0.9 condition, or implements a whole-event veto. Those requirements belong in the analyzer event-selection layer.

Direct unauthenticated retrieval of the authoritative [JME application documentation](https://cms-jme-jerc.docs.cern.ch/recommendations/jet-veto-maps/#application), [HIN jet ID](https://cms-hin-jet.docs.cern.ch/analysis/jetid), HIN `analysis/jetvetomap`, and [JMAR ID](https://cms-jme-jmar.docs.cern.ch/recommendations/jet-identification/#ak4ak8-chs) returned HTTP 401. Web-tool opens were also inaccessible; public GitHub repository searches did not locate their source projects. This is an access limitation, not evidence that a recommendation does not exist. The user-provided HIN prescription is the implementation requirement: Run 3 event veto, jet pT strictly above 15 GeV, applicable ID, CEF+NEF strictly below 0.9, any nonzero map value, and consistent data/MC application. Independent verification of the inaccessible page text remains outstanding.

The helper ID policies are:

| Detector eta interval | PbPb helper requested for OO | pp helper |
| --- | --- | --- |
| abs(eta)<=2.6 | NEF<0.9, MUF<0.8, CEF<0.8 | NHF<0.9, NEF<0.9, total constituents>1, MUF<0.8, CHF>0.01, charged multiplicity>0, CEF<0.8 |
| 2.6<abs(eta)<=2.7 | NEF<0.99, MUF<0.8, CEF<0.8 | NHF<0.9, NEF<0.99, MUF<0.8, charged multiplicity>0, CEF<0.8 |
| 2.7<abs(eta)<=3.0 | NEF<0.99 | NEF<0.99, neutral multiplicity>1 |
| 3.0<abs(eta)<=5.0 | NEF<0.9 | NHF>0.2, NEF<0.9, neutral multiplicity>10 |

Here charged multiplicity is CHM+CEM+MUM, neutral multiplicity is NHM+NEM, and total includes all five species. These are the exact category-count definitions in the selected HIN pp helper, so the implementation reproduces those sums. They should not be relabeled as a physically complete HF neutral count for an unrelated ID prescription. These helper policies are source evidence, not a claim that the initially supplied generic table has identical cuts.

## Immutable candidate identities

JECDatabase `master` was resolved once to commit `e825b8378a63df0cc4fb4c65dfe6088cbbb54a09`; all retrievals below use that immutable revision. The era/tag labels here are source directory names, not verified Git tags.

| Target | Source directory / era label | Exact ROOT filename | Size | SHA-256 |
| --- | --- | --- | ---: | --- |
| ppRef2024, pp RunJ | `jet_veto_maps/Summer24Prompt24_V1/` | `Summer24Prompt24_RunBCDEFGHI.root` | 484022 bytes | `ac1c38cdc92ab27f0362895591481c5e2e6a391dcb46359f385000d5105daedf` |
| OO2025 | `jet_veto_maps/Summer24Prompt25_V1/` | `Summer24Prompt25_RunCDEFG.root` | 481266 bytes | `a9b8f15da1bb417a784b647ac0199838642fba6ff6e324051b7309e759d68246` |

The exact source is [cms-jet/JECDatabase](https://github.com/cms-jet/JECDatabase): [pp file at the pinned commit](https://github.com/cms-jet/JECDatabase/blob/e825b8378a63df0cc4fb4c65dfe6088cbbb54a09/jet_veto_maps/Summer24Prompt24_V1/Summer24Prompt24_RunBCDEFGHI.root), [OO file at the pinned commit](https://github.com/cms-jet/JECDatabase/blob/e825b8378a63df0cc4fb4c65dfe6088cbbb54a09/jet_veto_maps/Summer24Prompt25_V1/Summer24Prompt25_RunCDEFG.root). Git blob IDs are respectively `b6a911ec640511254b08a254db51e01d90921457` and `1fd3b4a970074408483d0c1eb7f35cdc685eebff`.

The pp file was added by [6552021981eb7e0ec0bb25e87823c1e7c5d1338a](https://github.com/cms-jet/JECDatabase/commit/6552021981eb7e0ec0bb25e87823c1e7c5d1338a), dated 2025-07-17. The OO-nominated file was added by [17a7192932fd07daf7ad69920ab43f9610855bc7](https://github.com/cms-jet/JECDatabase/commit/17a7192932fd07daf7ad69920ab43f9610855bc7), dated 2026-07-01. Its [PR 226](https://github.com/cms-jet/JECDatabase/pull/226) describes the map as copied from `Winter25Prompt25_RunCDEFG.root`, within a campaign for 2025 data and Summer24 MC. That description does not explicitly establish OO reconstruction or special-run applicability.

## Run and production applicability

Neither ROOT file contains run-IOV metadata. The pp correctionlib companion also supplies no numeric run range. Thus the candidate numeric run ranges and full-sample coverage are **unknown**, not asserted from their filenames.

The authoritative public [CMS Tier-0 production record](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CompOpsTier0Policies) distinguishes the relevant special periods:

| Production-era record | Interval in that record | Bounded Goal 2 sample observed directly |
| --- | --- | --- |
| Combined Run2024J/HIRun2024A agent period | 387122–388260; preceding normal Run2024I period ends at 387121 | ppRef: 100 events, all run 387474 |
| Combined pORun2025/OORun2025/NeNeRun2025 agent period | 393610–394285, between normal Run2025C and Run2025D periods | OO: 1000 events, all run 394183 |

Those are combined production-agent intervals, not exact individual collision-system intervals or map validity intervals. The observed runs come from `goal2-native-hiforest-gic_ugf7/{ppRef_data_production,OO_data_production}.events.jsonl` in the retained 2026-10-28 Plan evidence. Their paths and SHA-256 are in the JSON. These bounded ledgers do not inventory the full campaigns. Normal-era labels `BCDEFGHI` and `CDEFG` therefore need explicit evidence extending the nominated maps to pp RunJ and OO; neither inclusion nor exclusion is proved solely by the labels.

The retained producer report pins `JesusVlz/cmssw@2400b134c3719accf688b3d936908c18343a1f15` and identifies pp CHS and OO constituent-subtracted PF targets. Exact producer-to-file attribution remains incomplete. An available pp data PSet retains jets only above 20 GeV; it is not established as the producer of the bounded older pp forest. Available MC PSets have raw and stored-pT retention thresholds of 15 GeV. These are relevant records, but they do not certify completeness for a veto defined with a potentially different pT convention or full detector eta coverage.

The intended derivation collection and radius are not encoded by either ROOT file and were not independently established for these productions. The proposed analyzer lookup uses detector E-scheme `jteta/jtphi`, while EEC analysis retains WTA geometry. The exact-production coordinate and raw/stored/new-corrected pT prescriptions still need applicable evidence; R=0.8 applicability is a separate requirement.

## Direct ROOT inspection

Both candidates were opened with PyROOT/ROOT 6.40.04. All objects are `TH2D`, cycle 1. Every histogram uses 82 variable eta bins spanning **[-5.191, 5.191]** on X and 72 uniform phi bins spanning **[-pi, pi] radians** on Y: 5904 regular bins. The OO titles are `#eta`/`#phi` except `jetvetomap_fpix`, whose axis titles are blank; pp axis titles are blank. The pp companion correctionlib names its inputs eta and phi. Full numeric edges and exact stored titles are retained in JSON.

| Histogram key | pp nonzero regular bins | OO nonzero regular bins |
| --- | ---: | ---: |
| `jetvetomap_hot` | 62 | 30 |
| `jetvetomap_cold` | 197 | 233 |
| `jetvetomap_hotandcold` | 259 | 263 |
| `jetvetomap_bpix` | 72 | 72 |
| `jetvetomap_fpix` | 8 | 26 |
| `jetvetomap_all` | 260 | 266 |
| `jetasymmetrymap` | 5904 | 5904 |
| `jetasymmetrymap_norm` | 5904 | 5904 |
| `jetpullsummap_nom` | 5904 | 5904 |
| `jetpullsummap_abs` | 5904 | 5904 |
| `jetpullsummap_ref` | 5904 | 5904 |
| `jetvetomap` | 164 | 169 |

All inspected regular contents are finite. All underflow/overflow contents are zero; these flow bins must not be interpreted as acceptance outside the physical axes. ROOT's upper phi edge is overflow, so periodic phi must be canonicalized before lookup. Cold-map veto bins contain -100, reinforcing the required `value != 0` predicate.

The user and pinned HIN helpers specify **`jetvetomap_all`**. However, the same pinned pp directory's [correctionlib companion](https://github.com/cms-jet/JECDatabase/blob/e825b8378a63df0cc4fb4c65dfe6088cbbb54a09/jet_veto_maps/Summer24Prompt24_V1/jetvetomaps.json.gz) recommends **`jetvetomap`** for data and MC. Its SHA-256 is `849cd1098814175ab58faac4834df0a353de63482f4017e22c0134ee7e840aec`. The two keys differ in 96 pp bins and 97 OO bins. They are not interchangeable. The requested key is recorded without silently substituting the correctionlib key; resolve the era-specific prescription before deployment.

The same detector-region map and eligible-jet event predicate are expected for data and corresponding MC under the supplied HIN recommendation, including signal and mixing populations. That does not imply equal event-loss fractions or shared event IDs. No event-loss, ID-efficiency, before/after data/MC, or mixing-closure measurement was run in this map-content audit.

## Reproduction and future deployment record

[inspect_veto_maps.py](inspect_veto_maps.py) reads local candidates, reports all keys/axes/content summaries, and can generate three-panel map comparisons. It performs no download or payload installation:

```bash
python3 validation/run3Jets/inspect_veto_maps.py \
  /tmp/run3-goal3-veto-audit/Summer24Prompt24_RunBCDEFGHI.root \
  /tmp/run3-goal3-veto-audit/Summer24Prompt25_RunCDEFG.root \
  --output /tmp/run3-goal3-veto-audit/direct_inspection.json \
  --plots-dir /tmp/run3-goal3-veto-audit/plots
```

The generated PNG paths and hashes are in the candidate JSON. These plots compare map contents only; they are not detector or event-selection validation. Copies are preserved under `validation/run3Jets/evidence/2026-10-28/` without copying a ROOT map.

The planned payload layout is:

```text
jetVetoMaps/
  OO2025/
    <approved-map>.root
    <approved-map>.metadata.json
  ppRef2024/
    <approved-map>.root
    <approved-map>.metadata.json
```

Every future deposited map must have a metadata record with source repository, immutable source commit/tag, SHA-256, era and numeric run coverage with full-sample coverage evidence, histogram key, jet collection, radius, coordinate convention, pT convention, data/MC applicability, validation plots, and approval evidence. Candidate records deliberately contain null/unestablished values rather than invented approvals. The assistant did not deposit any retrieved map. The user-uploaded copies have candidate metadata; no job package or production card was changed.

## Owner production clarification

After the audit, the owner reported that the ppRef data and MC ntuples were
produced with a 15 GeV threshold, then explicitly confirmed keeping 15 GeV. This owner-provided fact
is separate from the reference PSet with 20 GeV retained above. It does not
identify the archived PSet or resolve raw versus stored/corrected pT retention.
The retention threshold is now confirmed by the owner; the exact pT convention and other upstream cuts remain to verify. A 15 GeV bound
can suffice when it is in the veto's exact pT convention and covers all required
eta and eligible jets; a second upstream cut in another convention must also be
accounted for.
