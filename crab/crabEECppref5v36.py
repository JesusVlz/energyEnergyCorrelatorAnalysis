from WMCore.Configuration import Configuration
config = Configuration()

card='cardEECpp5v36TeV.input'
jobTag='eecAnalysis_ppref_2024-09-10'
inputList='ppRefForestFiles.txt'
outputFile=jobTag+'.root'
fileLocation='1'  # CERN

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+outputFile,'location='+fileLocation,'mixinglist=none']  # VERIFICAR punto 3
config.JobType.inputFiles = ['FrameworkJobReport.xml','eecOO.tar.gz',card]  # tarball propio, ver punto 1
config.JobType.outputFiles = [outputFile]
config.JobType.maxJobRuntimeMin = 120
config.JobType.maxMemoryMB = 1500

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'eecOOHistograms'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/jvelazqu/energycorr/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']