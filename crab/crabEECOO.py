import os
from pathlib import Path
from urllib.parse import urlsplit

from WMCore.Configuration import Configuration
config = Configuration()

baseDirectory = Path(__file__).resolve().parent
card = 'cardEECOO.input'
part = int(os.environ.get('OO_CRAB_PART', '1'))
chunk = int(os.environ.get('OO_CRAB_CHUNK', '1'))
if not 1 <= part <= 10 or not 1 <= chunk <= 2:
    raise ValueError('OO_CRAB_PART must be 1..10 and OO_CRAB_CHUNK must be 1..2')
splitTag = f'part{part:02d}_chunk{chunk:02d}'
jobTag = f'eecAnalysis_OO_mixedPerp_{splitTag}_Run3_v2'
inputList = baseDirectory / f'OOForestfiles/OOForestFiles_{splitTag}.txt'
mixingPool = baseDirectory.parent / 'mixingFileList/OO2025_MB_pool_pilot.txt'

# Vanderbilt is the default execution site. Add sites here to allow remote reads.
# Multiple whitelist entries are eligible sites, not ordered fallback choices.
additionalSites = ['T2_US_Purdue', 'T2_US_Nebraska']  # Optional: ['T2_US_Purdue', 'T2_US_Nebraska']

def readFileList(path):
    with path.open() as stream:
        return [line.strip() for line in stream if line.strip() and not line.lstrip().startswith('#')]

def logicalFileName(value):
    path = urlsplit(value).path if '://' in value else value
    offset = path.find('/store/')
    if offset < 0:
        raise ValueError('Expected a /store/ file name: ' + value)
    return path[offset:]

signalFiles = readFileList(inputList)
signalLFNs = {logicalFileName(value) for value in signalFiles}
poolLFNs = {logicalFileName(value) for value in readFileList(mixingPool)}
if not signalFiles or len(signalLFNs) != len(signalFiles):
    raise ValueError('Signal file list must be nonempty and contain no duplicates')
if len(signalFiles) > 10000:
    raise ValueError('One-file jobs require at most 10,000 signal files per CRAB task')
if len(poolLFNs) != 35:
    raise ValueError('The OO pilot expects 35 unique mixing-pool files')
if signalLFNs & poolLFNs:
    raise ValueError('Signal file list overlaps the reserved OO mixing pool')
outputFile=jobTag+'.root'
fileLocation='2'  # Vanderbilt

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = str(baseDirectory / 'PSet.py')
config.JobType.scriptExe = str(baseDirectory / 'compileAndRun.sh')
# The final analyzer argument is a legacy integer index; MixingFileList in the card selects the pool.
config.JobType.scriptArgs = ['card='+card, 'output='+outputFile, 'location='+fileLocation, 'mixinglist=0']
config.JobType.inputFiles = [str(baseDirectory / name) for name in ('FrameworkJobReport.xml', 'eecOO.tar.gz', card)]
config.JobType.outputFiles = [outputFile]
config.JobType.maxJobRuntimeMin = 120
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.userInputFiles = signalFiles
config.Data.splitting = 'FileBased'
# Isolate empty signal forests to their own jobs.
config.Data.unitsPerJob = 1
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'eecOOconesHistograms'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/jvelazqu/eec/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Vanderbilt'] + additionalSites
config.Site.storageSite = 'T2_CH_CERN'

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
