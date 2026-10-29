# Bounded Run-3 R4 calibration diagnostics

Preserve the inherited PbPb WTA definition: start reconstructed-jet correction
from rawpt, evaluate each level sequentially, use the selected WTA coordinates
for lookup, keep the WTA axis unchanged, append residuals only for data, then
apply analysis pT cuts. Root and CRAB cards differ; neither establishes the
published production configuration without retained production evidence.

`checkCalibrationInputs.cxx` checks the sequence and configured residual counts,
legacy lookup-coordinate preservation, disabled Run-3 JES/JER, and a future
Run-3 corrected-pT uncertainty input. That guard does not enable an uncertainty
payload or change nominal pT. Build it with current src objects and ROOT flags,
using `-fno-access-control` only for the regression fixture.

`compare_calibrations.py` reads bounded native R4 forests. Old and candidate
payloads share physical event IDs, reconstructed indices, native legacy matches,
weights, selection and WTA geometry, before analysis pT cuts. It saves a paired
JSONL ledger, ROOT histograms, pT/ratio/migration plots, and response mean/width
versus gen pT, eta and stored OO activity. Empty and low-statistics bins remain
explicit. Both energy powers use the same jet-level input sample; these are not
EEC-weighted response matrices or n=1/n=2 pair-closure results.

Candidate upstream L2 files are evaluated directly on rawpt as a **hypothesis**,
not installed as nominal. A separate coordinate-only diagnostic uses the same
installed payload and varies only lookup coordinates. V2's earlier CHS-labelled
provenance and the supplied “includes L1FastJet” statement do not by themselves
establish a standalone-L1 sequence. No L1 or rho/area input is invented. Inputs
with missing rawpt/WTA branches are recorded as unrun; stored jtpt or another
axis is not substituted. Data trigger/producer applicability must be established
separately. The data-only residual variation is isolated when compatible data
are actually available; residual changes do not change MC calibration.

Example, with pinned public upstream files retained outside the installed
payload directory:

```bash
python3 validation/run3Jets/compare_calibrations.py \
  --oo /path/to/actual_OO_MC.root --pp /path/to/actual_pp_MC.root \
  --upstream /path/to/pinned_HIN-Jets_snapshot \
  --output /tmp/new_r4_comparison --oo-events 10000 --pp-events 3000
python3 -m unittest discover -s validation/run3Jets -p 'test_*.py' -v
```

`veto_rule.py` tests the supplied **event** veto: any eligible reco jet with
pT>15 GeV, passing the applicable ID, CEF+NEF<0.9 and a nonzero map value rejects
the event. It requires explicit pT and detector-coordinate keys, a map lookup,
and an ID callback. It has no default map or ID. The standalone Python function is separate from the
explicit opt-in C++ analyzer integration described below. Era/run/collection/map-key applicability, threshold convention and
producer completeness must be established before signal/mixing activation.
Unavailable applicability means null event-loss fraction and after-veto maps
NOT_RUN, not zero loss. The plotted matched-jet before map is not a complete
eligible-veto-jet distribution. Native R8 remains an independent raw-only task.

The existing pair exporter selects stored jtpt and applies its documented
weights/selection. These diagnostics do not silently change or certify those
exports. Regenerate dependent exports/effective maps, held-out closures,
backgrounds, responses and downstream results when an approved calibration or
selection changes; track-only payloads require applicability checks rather than
automatic replacement. Historical artifacts retain their original meanings.

## Analyzer integration (Goal 3)

See [IMPLEMENTATION.md](IMPLEMENTATION.md) for the implemented, user-confirmed HIN
PbPb-for-OO and pp-for-ppRef ID, explicit event-veto configuration, conditional
QA and reproduction. [Candidate provenance](vetoMapSources.md) records the verified
user-uploaded ROOT files and unresolved production coverage. The earlier standalone
`veto_rule.py` remains a diagnostic; the analyzer now has the C++ event-level path.
