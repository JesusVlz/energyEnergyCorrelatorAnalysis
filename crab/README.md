# Running the code in crab

All the necessary configuration for running the code using the CRAB system is found from this folder. Below there is a simple example on how to run the analysis on CRAB and also detailed descriptions of all the files in this folder.

## Example commands

1. Create local tar ball:
    ```
    ./makeEECAnalysisTar.sh
    ```

2. Move the tar ball to lxplus:
   ```
   scp eec5TeV.tar.gz jviinika@lxplus.cern.ch:/afs/cern.ch/work/j/jviinika/private/eecAnalysis/CMSSW_12_4_8/src/eecAnalysiss
   ```

3. ssh to lxplus and set up the environment:
   ```
   ssh jviinika@lxplus.cern.ch
   cd /afs/cern.ch/work/j/jviinika/private/eecAnalysis/CMSSW_12_4_8/src/eecAnalysis
   cmsenv
   voms-proxy-init --voms cms
   ```

4. Modify the input card and crab configuration at lxplus:
   ```
   vim cardEECPbPb.input
   vim crabEECPbPb.py
   ```

5. Send the crab jobs:
   ```
   crab submit -c crabEECPbPb.py
   ```

## Executables

**compileAndRun.sh**

This is the script that is executed by CRAB. What is does is that it unzips the analysis tar ball, compiles the code, and runs the analysis.

## Crab configuration files

The configuration files given here are configured to sace the histograms to my (jviinika) EOS space. You will need to change that before you can use these. To do this, just put to variable config.Data.outLFNDirBase a path to a folder to which you have write right and possibly also change the config.Site.storageSite to a site to which to can write. Otherwise for each run you should set the variable infoString, which is the name of the run. It should describe the current configuration used. If you do not know where the input files for the analysis are located, setting fileLocation variable to '2' searches for the files. It is slightly slower to execute, but works with files located everywhere.

**crabEECPbPb.py** /
**crabEECPp.py**

Configuration file for running over the 2018 PbPb or 2017 pp dataset.

**crabEECPbPbMC.py** /
**crabEECPpMC.py**

Configuration file for running over the Pythia+Hydjet or Pythia8 samples.

## Analysis configuration files

These are the operational copies in `crab/`. Their paths are unchanged by
the Run-2 archive move. The Run-3 R=0.4 cards enable HIN Jet ID and the
whole-event veto as described in the [main README](../README.md). The CRAB configurations use local
card basenames, so resolve them from this directory or a staging directory
containing these copies. The archived root-level Run-2 local cards in
[`../cards/archive/run2/`](../cards/archive/run2/README.md) are independent
configurations; do not redirect the CRAB configurations to that archive.
See the [repository card inventory](../cards/README.md).

The configuration for the analysis is provided in card.input files. For example the collision system, data vs. MC running, and several analysis cuts can be defined in this file. The files contain comments on what each variable mean, so it should be straightforward to modify to your needs. The most important example files are.

**cardEECPbPb.input** /
**cardEECPp.input**

Configuration for running with 2018 PbPb or 2017 pp data.

**cardEECPbPbMC.input** /
**cardEECPpMC.input**

Configuration for running over the Pythia+Hydjet or Pythia8 samples.

## File lists

File lists are needed by the CRAB configuration. They contain names for all the files that are analyzed in the configured run. The following file lists are included:

**PbPbData2018_flowJetSkims_80and100triggers.txt**

These are the default files for the 2018 PbPb data. These are forested from /HIHardProbes/HIRun2018A-04Apr2019-v1/AOD primary dataset and either calo jet 80 or calo jet 100 triggers are required for each event in this skim.

**ppData2017_HighEGJet.txt**

Default file list for 2017 pp data. The files are forested from /HighEGJet/Run2017G-17Nov2017-v2/AOD primary dataset.

**newPbPbMcSkim.txt**

Default file list for Pythia+Hydjet simulation. The files are forested from /DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM dataset.

**ppMC2017_mergedFiles.txt**

Default file list for Pythia8 simulation. The files are forested from /QCD_pThat-15_Dijet_TuneCP5_5p02TeV_pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV_v1-v1/AODSIM dataset.

**triggerEfficiencyFiles2018.txt**

Small sample from /HIHardProbes/HIRun2018A-04Apr2019-v1/AOD dataset with no trigger selection and no tracks for trigger efficiency studies.

**PbPbData2018_MinBiasFiles.txt**

Small minimum bias sample for PbPb2018.

## Other files needed for Crab running

**PSet.py** /
**FrameworkJobReport.xml**

Dummy files. Not used, but need to be there to keep CRAB happy.
