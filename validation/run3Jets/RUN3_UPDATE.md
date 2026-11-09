# Run-3 jet-pT calibration, Jet ID and QA update (8 November 2026)

This review branch uses the owner-supplied HIN miniPOG Jet ID table and the current
Run-3 R=0.4 jet-energy payloads relevant to the EEC analysis.

## Jet-energy corrections

The analysis still treats JEC as a correction of reconstructed jet pT. For the
2024 ppRef L1FastJet level, the calibration formula additionally requires the
event pileup density and measured jet area internally. Following the calibration
reference implementation, these are read from:

- `ggHiNtuplizer/EventTree/rho`
- `ak4PFJetAnalyzer/t/jtarea`

Detector `jteta` selects the calibration row. WTA eta/phi remain the EEC
analysis axis.

The active chains are:

| sample | active chain |
| --- | --- |
| OO 2025 data | rawPt -> Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt |
| OO 2025 MC | rawPt -> Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt |
| pp 2024 MC | rawPt -> Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt -> Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt |
| pp 2024 data | rawPt -> Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt -> Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt -> Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt |

The pp ordering follows the 2024 ppRef calibration implementation in
`thenicholasbarnett/HIN-L2Residuals-JER`: L1FastJet is applied before
L2Relative to both data and MC, while the residual is appended only for data.

No OO residual is applied in the adopted prescription.

Payload layout:
- `jetEnergyCorrections/OO2025/`
- `jetEnergyCorrections/pp2024/`

The files are pinned to `thenicholasbarnett/HIN-Jets@c64293a132af388af1827044779af36e57f884d2`.

## Jet ID and trackMax

`Run3JetID 1` supplements the existing strict
`MinMaxTrackPtFraction < trackMax/rawPt < MaxMaxTrackPtFraction` requirement.
The nominal bounds remain 0.01 and 0.98. The denominator is raw jet pT.

OO follows the supplied PbPb Run-3 ID table. pp follows the supplied Run-3
TightLeptonVeto table. Generator jets do not receive reconstructed PF ID or
trackMax requirements.

The whole-event detector veto remains a separate selection and does not gain a
trackMax requirement.

## QA output

`Run3JetQA 1` stores:
- baseline / trackMax / Jet ID / combined / vetoed selection-stage jet-pT spectra,
- `trackMax/rawPt`,
- fixed-final-jet calibration spectra at raw, L1FastJet, L2Relative and
  (pp data only) L2Residual stages,
- a snapshot of the configured jet-veto map.

Example:

```bash
root -l -b -q 'plotting/plotRun3JetQA.C("OO_output.root","oo50to80",50,80)'
root -l -b -q 'plotting/plotRun3JetQA.C("pp_output.root","pp",-1,0)'
root -l -b -q 'plotting/plotRun3VetoMap.C("jetVetoMaps/OO2025/Summer24Prompt25_RunCDEFG.root","OO_veto.png")'
root -l -b -q 'plotting/plotRun3VetoMap.C("jetVetoMaps/pp2024/Summer24Prompt24_RunBCDEFGHI.root","pp_veto.png")'
```

Existing outputs made before this update cannot recover the pre-cut QA stages and
should not be mixed with responses produced using the new additive jet selection.
