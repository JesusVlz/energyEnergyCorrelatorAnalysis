from WMCore.Configuration import Configuration
from datetime import date
config = Configuration()

iPart = '2'
system = 'RecoGen'
comment = 'eWSquared_nomSmear_jetPtW_recoRef'
today = date.today().strftime('%Y-%m-%d')
infoString = 'Pythia8_pfJets_wtaAxis_' + comment + '_part' + iPart + '_' + today

card='responseMatrixVariationCardPpMC.input'
inputFile='ppMC2017_newForest_part' + iPart + '.txt'
output='ppMC2017_' + system + '_' + infoString + '.root'
fileLocation='1'  # Locations: 0 = Purdue, 1 = CERN, 2 = Vanderbilt,  3 = Search with xrootd

config.section_("General")
config.General.requestName = 'ppMC2017_' + system + '_' + infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','eec5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 120
config.JobType.maxMemoryMB = 1200

config.section_("Data")
config.Data.userInputFiles = open(inputFile).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'eecPpMCHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

