# Run-3 jet-pT calibration, Jet ID and QA update (8 November 2026)

This review branch uses the owner-supplied Run-3 Jet ID table, veto maps, and
current R=0.4 JEC payloads relevant to the EEC analysis.

## Jet-energy corrections

The active ppRef production cards intentionally set:

```text
Run3JEC 1
Run3UsePpL1FastJet 0
```

The current EEC ppRef forests do not save the event rho input required by the
L1FastJet payload. Therefore L1FastJet is available in the framework but is not
applied to current production. No rho=0 or other guessed value is used.

Active production chains:

| sample | active chain |
| --- | --- |
| OO 2025 data/MC | rawPt -> Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt |
| pp 2024 MC | rawPt -> Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt |
| pp 2024 data | rawPt -> Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt -> Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt |

If a future ppRef forest contains the required inputs, setting
`Run3UsePpL1FastJet 1` activates:

- `ggHiNtuplizer/EventTree/rho`
- `ak4PFJetAnalyzer/t/jtarea`
- `rawPt -> L1FastJet -> L2Relative -> [data residual]`

The switch is rejected outside ppRef2024 or when `Run3JEC 0`.

Detector `jteta` selects the calibration row. WTA eta/phi remain the EEC
analysis axis. The L1 payload is retained under `jetEnergyCorrections/pp2024/`
for future use.

## Jet ID and trackMax

`Run3JetID 1` supplements the existing strict
`MinMaxTrackPtFraction < trackMax/rawPt < MaxMaxTrackPtFraction` requirement.
The nominal bounds remain 0.01 and 0.98, with raw jet pT in the denominator.

OO follows the supplied PbPb Run-3 ID table. pp follows the supplied Run-3
TightLeptonVeto table. The whole-event detector veto remains separate and does
not gain a trackMax requirement.

## QA output

`Run3JetQA 1` stores:
- baseline / trackMax / Jet ID / combined / vetoed jet-pT spectra,
- `trackMax/rawPt`,
- fixed-final-jet calibration spectra,
- a snapshot of the configured veto map.

The calibration histogram always has labeled slots for raw, L1FastJet,
L2Relative and data residual. With `Run3UsePpL1FastJet 0`, the L1 slot is
intentionally empty and the plotting macro automatically omits it.

Existing outputs made before this update should not be mixed with responses
produced using the updated calibration/selection configuration.
