# Code to construct correlations

This folder has the code that constructs correlations from forests. The code in compiled with the Makefile in the parent directory. Below is given a brief explanation what each class does.

**EECAnalyzer.cxx**

The main correlation code can be find from here. The code loops over a list of files given into it. It goes through all the events in the files, and analyzes them. It fills track, inclusive jet and energy-energy correlator histograms. Most of the configuration for the analysis can be read from an configuration card.

**EECHistograms.cxx**

All the histograms that are filled by EECAnalyzer are managed by this class. If you ever need to change binning, add histograms or anything else regarding the histograms, it can be done in this class. THnSparses are used for most of the analysis histograms, since they can conveniently contain several binning options.

**ForestReader.cxx**

A class template defining what kind of functionality a class reading the trees from a forest must implement. There is a separate class reading the trees since the analysis should not depend on the exact implementation of the forest. If the forest changes, one should write a new class inheriting from the ForestReader class reading it and no changes should be made in the main analysis code.

**HighForestReader.cxx**

A forest reader that can read the standard CMS heavy ion forests.

**GeneratorLevelForestReader.cxx**

A forest reader that reads generator level jet and track information instead of reconstructed values.

**ConfigurationCard.cxx**

A class that extracts information from specifically formatted text file (.input) and stores it. This allows to change input parameters for the analysis without touching the source code. It also creates a reference on the parameter set that is used for each run that is stored in the root file with the correlations.

**JetCorrector.cxx**

Class for doing jet energy correction on the fly. The code is originally from https://twiki.cern.ch/twiki/pub/CMS/HiJetReco2019/JetCorrector.h. It has been separated to cxx and h files and some debugging options have been added.

**JetUncertainty.cxx**

Class for estimating an uncertainty for the jet energy. The code is originally from https://twiki.cern.ch/twiki/pub/CMS/HiJetReco2019/JetUncertainty.h. It has been separated to cxx and h files and some debugging options have been added.

**TrackingEfficiencyInterface.cxx**

Class interface that defines what functions tracking efficiency classes for new data must implement. Used to seemslessly switch between pp and PbPb tracking corrections in DijetAnalyzer.

**trackingEfficiency2018PbPb.cxx**

Class managing tracking efficiency correction for 2018 PbPb data. The code is originally from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HITracking2018PbPb. It has been slightly modified for this analysis.

**trackingEfficiency2017pp.cxx**

Class managing tracking efficiency correction for 2017 pp data. The code is originally from https://twiki.cern.ch/twiki/bin/view/CMS/HiTracking2017pp. It has been slightly modified for this analysis.

## Slim OO/ppRef forests

For DataTypes 9-12, the five PbPb-only MiniAOD track fields
`trkNormChi2`, `trkNLayers`, `trkNHits`, `pfEcal`, and `pfHcal`
are optional in HighForestReader. Present fields are still read; absent
fields are not bound, and their getters report an error if called.
The OO/ppRef track cuts return before those PbPb-only quality cuts.
PbPb's existing reader requirements are preserved. UnfoldingForestReader
already treats these fields as optional.

The slim track contract retains `nTrk`, `trkPt`, `trkEta`, `trkPhi`,
`trkPtError`, `highPurity`, the four first-vertex impact-parameter fields,
and `trkCharge` (for response studies). Native MiniAOD vector types and
event-tree entry ordering must be preserved.

This compatibility change does not relax the current native jet input
contract: `rawpt`, `trackMax`, measured `jtarea` when L1 is enabled,
and the ten `jtPf*` composition branches must be genuinely produced.
Old `chargedMax`/`jtptUncorrected` schemas are not automatically aliased.

For ppRef R4, optional L1FastJet uses scalar
`ggHiNtuplizer/EventTree/rho` from `fixedGridRhoFastjetAll` and measured
`ak4PFJetAnalyzer/t/jtarea`. Enable it with `Run3JEC 1` and
`Run3UsePpL1FastJet 1` only on forests containing valid, aligned inputs.
The OO eta-binned rho profile is not this scalar, and current Run-3
calibration/QA remains restricted to native R4.

The matching producer configuration and a checker that enumerates all
ROOT trees/branches are maintained on
[the CMSSW EEC review branch](https://github.com/JesusVlz/cmssw/tree/eec-ppref-l1fastjet-20261109/HeavyIonsAnalysis/Configuration).

For ppRef MiniPOG Jet ID, new forests also save Int_t[nref] `jtPfCM`
and `jtPfNM` from the aggregate PAT charged/neutral multiplicities.
Their sum supplies the constituent count, including HF constituents.
The five category sums remain separate for the pinned HINpp prescription.
Legacy forests without the pair retain the unresolved-lower-bound error;
reforest to supply these values. Partial pairs and schema changes across
a chain are rejected. The producer update adds two integers per ppRef jet.

A ROOT-based regression test covers legacy unresolved/resolved counts,
HF-inclusive exact counts, genuine one-constituent rejection, invalid pairs,
overflow, and chain schema consistency:

```bash
g++ -std=c++17 $(root-config --cflags) validation/test_full_pf_multiplicities.cxx src/Run3JetSelection.cxx $(root-config --libs) -o /tmp/test_full_pf_multiplicities
/tmp/test_full_pf_multiplicities
```
