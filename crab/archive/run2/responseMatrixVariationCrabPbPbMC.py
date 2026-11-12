from WMCore.Configuration import Configuration
from datetime import date
config = Configuration()

iPart = '2'
system = 'RecoGen'
centralityShift = '4pC'
comment = 'eWSquared_nomSmear_jetPtW_recoRef'
today = date.today().strftime('%Y-%m-%d')
infoString = 'akFlowJets_' + centralityShift + '_cutBadPhi_' + comment + '_part' + iPart + '_' + today

card='responseMatrixVariationCardPbPbMC.input'
output='PbPbMC2018_' + system + '_' + infoString + '.root'
inputFile='pythiaHydjet2018_miniAODforest_part' + iPart + '.txt'
fileLocation='2'  # Locations: 0 = Purdue, 1 = CERN, 2 = Vanderbilt,  3 = Search with xrootd

config.section_("General")
config.General.requestName = 'PbPbMC2018_' + system + '_' + infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','eec5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 400
config.JobType.maxMemoryMB = 1600

config.section_("Data")
config.Data.userInputFiles = open(inputFile).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'eecPbPbMCHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

