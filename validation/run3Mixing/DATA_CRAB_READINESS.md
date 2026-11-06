# Run-3 R4 data mixing / CRAB readiness

Status: **not ready for submission**. This audit covers OO split and inclusive
data and ppRef 5.36 TeV data. No data task was submitted by this audit.

## Selected matching policy

The owner's latest instruction selects `EventMatchForMixing 0` for every channel:
Vz plus hiBin. This supersedes the earlier request for Vz-only inclusive/pp
matching. Inclusive OO still has no centrality selection cut. Matching two
events in hiBin does not define a cut on the signal sample's centrality.

For ppRef, `HighForestReader` does not bind the stored hiBin branch: it sets
`fHiBin=-1`, and `GetHiBin()` returns the constant 1 for matching. Therefore
option 0 already acts as Vz-only matching for ppRef. This is a reader convention,
not a claim that every upstream pp forest stores an identical hiBin value.

The six relevant R4 data cards now use option 0: `cardEECOO.input`,
`cardEECOOInclusiveR04.input`, `cardEECpp5v36TeV.input`,
`cardEECpp5v36TeVR04.input`, and the OO/ppRef copies under `crab/`.
Only the matching option/comments changed. Background enablement, thresholds,
and calibration settings still need a dedicated production configuration.
The exact card hashes and numeric-setting check are saved in
`results/2026-11-03-data-crab-readiness/option0_card_check.json`.

Existing options in `src/EECAnalyzer.h` are:

| Option | Variables | Initial unbinned tolerances |
| --- | --- | --- |
| 0 | Vz and hiBin | 0.5 cm; identical hiBin |
| 1 | Vz and HF energy | 0.5 cm; 1% above 100, otherwise 1 stored HF unit |
| 2 | Vz and background multiplicity | 0.5 cm; 1% above 100, otherwise 1 particle |

`BinnedMixing` selects numerical windows versus coordinate bins; it is not an
additional activity variable. The current general matchers automatically widen
tolerances after full pool scans and have no fixed search cap. The bounded OO
GenGen pilot is a separate path and does not establish data readiness.

## What the analysis note actually specifies

Reference: `/eos/user/j/jvelazqu/HIN/Plan/AN2023_081_v13.pdf`.
Page numbers below are the printed AN page numbers.

- **Section 7.4, p. 41, lines 618–622:** minimum-bias donors; vertex separation
  within 0.5 cm and identical hiBin (a 0.5%-wide centrality bin); two donors per
  signal event; each donor cone at the signal jet's eta and phi. This is the
  PbPb data prescription corresponding to option 0's initial window.
- **Section 7.4.1, pp. 42–43, lines 656–675:** embedding the hard event biases
  hiBin matching toward more active minimum-bias backgrounds. The generator
  study therefore matches HYDJET charged-particle multiplicity, pT>0.7 GeV and
  |eta|<2.4, within 1% (within one particle below 100), retaining the 0.5 cm Vz
  window. This generator-only definition is not a data observable.
- **Section 7.4.2, p. 49, lines 684–708:** adequate donor statistics matter;
  repeated use of the same events can amplify fluctuations. It does not justify
  arbitrarily widening the matching windows to force a match.
- **Section 9.5, p. 255, lines 1404–1407:** nominal pp results have no background
  subtraction; reflected-cone subtraction supplies a systematic comparison.
  Thus ppRef option-0 mixing is a requested Run-3 study, not an AN pp mandate.
- **Section 4, p. 10, lines 223–226:** the statement that pp needs only Vz
  matching concerns MC-to-data distribution reweighting, not donor matching.

Applying the PbPb prescription to OO requires Run-3 support/closure checks. The
AN's fixed matching window does not validate the code's automatic widening.
No new numerical pp tolerance has been selected by this audit.

## Remaining launch blockers and next actions

1. **Independent donor forests:** all 35 listed OO MB pilot donor LFNs occur in
   the 132,898-file signal list. Reserve donors outside the campaign's signal
   list (or verify a distinct pool), and reject overlaps after canonicalizing
   local/XRootD aliases. No ppRef2024 donor list is currently configured.
   Also exclude duplicate physical run/lumi/64-bit event identities.
2. **General data matcher:** enforce forest independence, preserve full event
   identities, use the total donor-chain count, and bound searches with explicit
   unmatched-event accounting. The binned hiBin vector path also needs repair
   before it can be used. Keep unbinned option 0 for the planned cards.
3. **Production cards:** existing CRAB cards use `BackgroundMethods 3` and
   `mixinglist=none`; they run reflected/perpendicular cones, not mixed cones.
   Bind explicit validated pools and a mixed/perpendicular configuration, with
   the agreed R4 bins, track threshold and both powers, before a pilot.
   The pp cards also retain `MegaSkimMixing 1`; use the regular-forest mode for
   ordinary donor forests and verify the pp ZeroBias donor selection. The CRAB
   baseline copies lack the dedicated cards' strict `MinPairTrackPt 2` setting.
4. **Frozen worker package:** `crab/eecOO.tar.gz` is stale and lacks the current
   Run-3 jet-selection helper and veto maps. Rebuild from verified source and
   payloads; propagate compile/analyzer errors and validate output plus ledgers.
5. **Data acceptance and runtime:** establish certified-lumi/skim provenance,
   trigger/retention applicability, valid submission credentials and pool
   staging/read access. Measure pp support and CPU in a bounded real-data pilot
   before deciding whether a looser window is needed.
6. **Pilot gate:** verify forest/physical-event exclusions, matching distances,
   two distinct donors, matched/unmatched exposure, source-event correction
   coordinates and reuse, plus runtime. Keep full calibration and pair-payload
   qualification separate from successful job execution.

Detailed audit evidence is under
`results/2026-11-03-data-crab-readiness/`: `analyzer_audit.md`,
`input_card_audit.md`, `crab_audit.md`, and their machine-readable records.

Verification: all six cards pass the numeric-setting comparison; only the pp
matching option changed numerically. Existing OO pilot matching and two-file
chain regressions pass under ROOT 6.32.13. Those regressions validate reusable
MC pilot mechanisms, not the ordinary data path or CRAB execution. In particular,
the chain fixture still demonstrates legacy count 2 versus total-chain count 5.
