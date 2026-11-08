# Run-3 calibration, Jet ID and stage QA update (8 November 2026)

This review branch follows the calibration table and HIN miniPOG Jet ID screenshot
supplied by the owner. The live CERN Jet ID page was inaccessible to this session;
`MiniPOG20261108` records the date supplied, NOT a claimed website release date.
Reference: https://cms-hin-jet.docs.cern.ch/analysis/jetid .
The older GitHub helpers remain available as a separately named cut table.

## Payloads and activation

All four files in `jetEnergyCorrections/run3/HIN-Jets_c64293a/` are pinned to
`thenicholasbarnett/HIN-Jets@c64293a132af388af1827044779af36e57f884d2`.
Verify them with `python3 validation/run3Jets/sync_review_payloads.py`.
The script checks Git blob hashes and record sizes, and reports SHA-256.

| Sample | `Run3JEC 1` actual chain |
| --- | --- |
| OO data and MC | `Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt` |
| ppRef MC | `Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt` then `Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt` |
| ppRef data | the MC chain then `Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt` |

There is no OO residual in this prescription. Do not substitute the upstream
centrality-tagged 80-90% file. Residuals are never applied to MC.
All levels start from native rawpt and use the preceding level's momentum.
`Run3JEC 1` uses detector jteta/jtphi while EEC geometry stays WTA.
Main analysis, reflected-jet QA, responses, closure and recomputed veto momentum
use the same calibration helper. Generator jets receive no reconstructed JEC.

The OO R4 cards enable the current payload. **ppRef cards retain `Run3JEC 0`
explicitly: the applicable rho source for the actual production has not been
established.** This is not a completed pp production calibration deployment.
To activate the implemented complete pp chain, set `Run3JEC 1` and supply:

```
Run3RhoTree <verified-tree-path>
Run3RhoBranch <verified-event-rho-scalar-branch>
Run3JetAreaBranch <verified-native-per-jet-area-branch>
```

These placeholders are NOT runnable branch names. Inspect the actual forest and
producer before filling them. The reader accepts Float_t/Double_t scalar rho and
Float_t/Double_t area[nref], checks types, count/entry alignment and chain
transitions, and rejects missing/nonfinite inputs or nonpositive area. It never
uses zero rho, a geometric area, or stored jtpt as a replacement for missing
calibration inputs. Entry alignment is checked; the chosen rho definition,
units, producer provenance and physics applicability still require confirmation.
A vector/differential rho producer needs a separately specified prescription;
this code deliberately does not invent how to reduce it to a scalar.

## Jet ID and the independent trackMax selection

`Run3JetID 1` adds ID to the unchanged strict requirement
`MinMaxTrackPtFraction < trackMax/rawpt < MaxMaxTrackPtFraction`.
The configured nominal bounds are 0.01 and 0.98. The denominator is NOT corrected
pT. Both cuts apply to reconstructed analysis, reflected-QA and response jets.
The event-veto eligibility scan does NOT gain a trackMax requirement.
Generator jets do not receive either reconstructed quality cut.

`Run3JetIDTable MiniPOG20261108` is the new default. OO follows the PbPb table:
NEF < 0.90/0.99/0.99/0.40 in |eta| <=2.6 / 2.6-2.7 / 2.7-3 / 3-5;
MUF<0.80 and CEF<0.80 in |eta|<=2.7, with no hadronic/multiplicity cuts.
pp follows TightLeptonVeto: NHF<0.99 centrally, NHF<0.90 at 2.6-2.7,
NHF<0.99 at 2.7-3, no NHF requirement beyond 3; NEF<0.90/0.99/0.99/0.40;
central CHF>0.01, constituents>1, charged multiplicity>0 through 2.7;
neutral multiplicity>1 at 2.7-3 and >10 beyond 3; MUF/CEF<0.80 through 2.7.
`HINc64293a` retains the historical helper table for explicit comparisons only;
selecting it does not revert the additive trackMax policy.

Native five-category neutral/total counts can omit HF constituents. For the
new generic pp table, a lower bound above a threshold establishes a pass;
an unresolved threshold raises a visible error. Full forward multiplicity
support must be confirmed; do not silently replace the new table by the old one.
This matters to the event veto even when EEC jets are central.

## Diagnostic output

`Run3JetQA 1` saves these objects in the ordinary `./eecAnalysis` ROOT output:

- `run3JetQA/selectionPtWeighted` and `selectionPtCounts`: common final calibrated
  pT, six quality categories, centrality. Categories include baseline, trackMax
  only, ID only, their intersection, ID+event-veto, and intersection+event-veto.
  Filled after trigger/standard event cuts, before event-veto rejection. The
  common EEC-axis eta, calibrated pT and matching requirements remain in place.
- `run3JetQA/calibrationPtWeighted`: raw, after L1, after L2, after data residual,
  for exactly the same final-selected reconstructed jets. Absent/inactive levels
  remain empty. This is distinct from a quality-cut comparison.
- `run3JetQA/trackMaxOverRawPt`: the actual ratio versus RAW pT and centrality.
- `run3JetQA/configuredJetVetoMap`: a detached copy of the exact configured mask,
  plus file/key and selection/calibration metadata. This is not occupancy.

These diagnostics are for reconstructed analysis-jet modes; GenReco/GenGen
outputs have no reconstructed analysis-jet stage entries. pp uses centrality
-0.5. pT >1000 GeV goes to overflow. QA does not change analysis selection or
weights. Stored Sumw2 is a sum of per-jet squared weights, not event-cluster
covariance. Ratio plots intentionally show central values without an unsupported
independent or binomial error model. Use event-level resampling for uncertainty.

From the repository root, after producing a new ROOT output:

```
root -l -b -q 'plotting/plotRun3JetQA.C("output.root","oo50to80",50,80)'
root -l -b -q 'plotting/plotRun3JetQA.C("pp_output.root","pp",-1,0)'
root -l -b -q 'plotting/plotRun3VetoMap.C("jetVetoMaps/OO2025/Summer24Prompt25_RunCDEFG.root","OO_veto.png")'
root -l -b -q 'plotting/plotRun3VetoMap.C("jetVetoMaps/ppRef2024/Summer24Prompt24_RunBCDEFGHI.root","pp_veto.png")'
```

Cumulative plot: baseline -> trackMax -> trackMax+ID -> trackMax+ID+event veto.
Independent plot reproduces the screenshot's comparison: baseline, trackMax
only, ID only, ID+event veto. Do not normalize each curve to its own jet/event
count: that would conceal rejection losses. Both plot families use one exposure.
Existing final-only outputs cannot recover pre-cut or pre-veto spectra; rerun.

Map plots preserve signed bin content and detector eta/phi axes. Plot the original
payload or one single-job snapshot; `hadd` adds TH2 snapshots, so a merged mask
snapshot does NOT retain original values. Ordinary QA spectra may be summed.
The existing `inspect_veto_maps.py` remains available for all-key comparisons.
Map IOV/production-coverage status is unchanged by adding a plotting macro.

Rebuild and rerun responses/corrections made with the previous ID-only selection;
structural binning compatibility alone does not establish selection compatibility.
No production-forest replay or physics closure is claimed by synthetic tests.
