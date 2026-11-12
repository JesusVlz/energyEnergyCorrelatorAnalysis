"""Shared, checked configuration for the four Run-3 EEC CRAB campaigns."""

import os
import re
from pathlib import Path
from urllib.parse import urlsplit


CRAB_DIR = Path(__file__).resolve().parent
REPO = CRAB_DIR.parent
EOS_PREFIXES = {
    'oo-mc': 'root://eosuser.cern.ch//eos/user/j/jvelazqu/Foresting/OOMC/',
    'ppref-mc': 'root://eosuser.cern.ch//eos/user/j/jvelazqu/Foresting/ppMC/',
}
SAMPLES = {
    'oo-data': ('oo', 9, 'OO_Data_split_sameEventCones.input', 'eecOOData'),
    'ppref-data': ('ppref', 11, 'ppRef_Data_inclusive_sameEventCones.input', 'eecPpRefData'),
    'oo-mc': ('oo', 10, 'OO_MC_split_none.input', 'eecOOMC'),
    'ppref-mc': ('ppref', 12, 'ppRef_MC_inclusive_none.input', 'eecPpRefMC'),
}


def read_list(path):
    if not path.is_file():
        raise ValueError(f'Missing input list {path}; see crab/README.md')
    lines = [s.strip() for s in path.read_text().splitlines()
             if s.strip() and not s.lstrip().startswith('#')]
    if not lines or len(lines) != len(set(lines)):
        raise ValueError(f'Empty or duplicate input list: {path}')
    return lines


def lfn(value):
    path = urlsplit(value).path if '://' in value else value
    offset = path.find('/store/')
    if offset < 0:
        raise ValueError(f'Expected /store/ in {value}')
    return path[offset:]


def input_list(sample, part=1, chunk=1):
    override = os.environ.get('EEC_INPUT_LIST')
    if override:
        selected = Path(override).resolve()
        if not selected.is_relative_to(CRAB_DIR / 'inputLists'):
            raise ValueError('EEC_INPUT_LIST must be a generated list under crab/inputLists/')
        return selected
    if sample == 'oo-data':
        if not (1 <= part <= 10 and 1 <= chunk <= 2):
            raise ValueError('OO data part=1..10, chunk=1..2')
        return CRAB_DIR / f'OOForestfiles/OOForestFiles_part{part:02d}_chunk{chunk:02d}.txt'
    if sample == 'ppref-data':
        if not (1 <= part <= 4 and chunk == 1):
            raise ValueError('ppRef data part=1..4, chunk=1')
        return CRAB_DIR / f'ppRefForestfiles/ppRefForestFiles_part{part:02d}_chunk01.txt'
    if part != 1 or chunk != 1:
        raise ValueError('MC has a single input list; part and chunk must be 1')
    return CRAB_DIR / 'inputLists' / ('OOMC_EOS.txt' if sample == 'oo-mc' else 'ppRefMC_EOS.txt')


def card_for(sample, choice=None):
    system, dtype, default, _ = SAMPLES[sample]
    if choice is None:
        choice = f'cards/R4/measurement/{default}'
    path = (REPO / choice).resolve() if not Path(choice).is_absolute() else Path(choice).resolve()
    if not path.is_relative_to(REPO / 'cards' / 'R4') or path.suffix != '.input':
        raise ValueError('Select a runnable .input card under cards/R4 (not a template)')
    if not path.is_file():
        raise ValueError(f'Missing card: {path}')
    content = path.read_text()
    if re.search(r'@[A-Z][A-Z0-9_]*@', content):
        raise ValueError(f'Unresolved template token in {path}')
    values = {}
    for line in content.splitlines():
        fields = line.split('#', 1)[0].split()
        if len(fields) >= 2:
            values[fields[0]] = fields[1]
    if int(values.get('DataType', -1)) != dtype:
        raise ValueError(f'Card DataType differs from {sample}: {path}')
    if sample.endswith('mc') and int(values.get('BackgroundMethods', '0')) & 4:
        raise ValueError(
            'Mixed-cone MC is not supported by this CRAB campaign. The OO HIJING '
            'pool A list needs event-level donor reservations and ledger outputs; '
            'use the existing strict OO GenGen closure workflow in condor/README.md.')
    if sample == 'ppref-data' and values.get('Run3UsePpL1FastJet', '0') != '0':
        raise ValueError('Vanderbilt ppRef forests have not been certified for L1FastJet inputs')
    return path, values


def task_name(sample, card, part, chunk, tag):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,23}', tag):
        raise ValueError('EEC_RUN_TAG must be 1..24 letters, digits or underscores, starting with a letter')
    stem = re.sub('[^A-Za-z0-9_]', '_', card.stem)
    label = SAMPLES[sample][3]
    name = f'{label}_{stem}_{tag}'
    if sample.endswith('data'):
        name += f'_p{part:02d}c{chunk:02d}'
    if len(name) > 100:
        raise ValueError('CRAB request name exceeds 100 characters: ' + name)
    return name


def files_for(sample, part, chunk, values):
    lines = read_list(input_list(sample, part, chunk))
    if sample.endswith('data'):
        if len(lines) > 10000:
            raise ValueError('One-file data jobs must fit the 10,000-job task limit')
        normalized = [lfn(s) for s in lines]
        if len(set(normalized)) != len(lines):
            raise ValueError('Duplicate signal LFNs')
        for entry in lines:
            if '://' in entry and urlsplit(entry).scheme != 'root':
                raise ValueError('Expected a /store LFN or XRootD PFN: ' + entry)
        # CRAB accepts LFNs. The analyzer's location=2 then consistently
        # prefixes Vanderbilt's XRootD endpoint for the actual read, even
        # when an older versioned list contains a different redirector PFN.
        lines = normalized
        if int(values.get('BackgroundMethods', '0')) & 4:
            selected = values.get('MixingFileList')
            selected_path = (REPO / selected).resolve() if selected else None
            if (not selected_path or not selected_path.is_relative_to(REPO / 'mixingFileList')
                    or not selected_path.is_file()):
                raise ValueError('Mixed-cone card requires a real packaged mixing list')
            mixed_donors = [lfn(s) for s in read_list(selected_path)]
            if len(set(mixed_donors)) != len(mixed_donors) or set(normalized) & set(mixed_donors):
                raise ValueError('Selected mixing list has duplicates or overlaps signal files')
    else:
        prefix = EOS_PREFIXES[sample]
        if any(not s.startswith(prefix) or not s.endswith('.root') for s in lines):
            raise ValueError(f'{sample} list contains non-EOS or unexpected files')
        if (len(lines) + 1) // 2 > 10000:
            raise ValueError('MC input list exceeds 10,000 jobs with two files/job')
    return lines


def create_config(sample):
    from WMCore.Configuration import Configuration

    part = int(os.environ.get('OO_CRAB_PART' if sample == 'oo-data' else 'PPREF_CRAB_PART', '1'))
    chunk = int(os.environ.get('OO_CRAB_CHUNK' if sample == 'oo-data' else 'PPREF_CRAB_CHUNK', '1'))
    if sample.endswith('mc'):
        part = chunk = 1
    card, values = card_for(sample, os.environ.get('EEC_CARD'))
    files = files_for(sample, part, chunk, values)
    tag = os.environ.get('EEC_RUN_TAG', 'R4v1')
    name = task_name(sample, card, part, chunk, tag)
    output = name + '.root'
    config = Configuration()
    config.section_('General')
    config.General.requestName = name
    config.General.workArea = str(CRAB_DIR / 'crab_projects' / 'run3')
    config.General.transferOutputs = True
    config.General.transferLogs = True
    config.section_('JobType')
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = str(CRAB_DIR / 'PSet.py')
    config.JobType.scriptExe = str(CRAB_DIR / 'compileAndRun.sh')
    config.JobType.scriptArgs = [f'card={card.name}', f'output={output}',
                                 'location=' + ('2' if sample.endswith('data') else '1'), 'mixinglist=0']
    config.JobType.inputFiles = [str(CRAB_DIR / 'FrameworkJobReport.xml'),
                                 str(CRAB_DIR / 'eecOO.tar.gz'), str(card)]
    config.JobType.outputFiles = [output]
    config.JobType.maxJobRuntimeMin = 120
    config.JobType.maxMemoryMB = 2500
    config.section_('Data')
    config.Data.userInputFiles = files
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1 if sample.endswith('data') else 2
    config.Data.totalUnits = len(files)
    config.Data.outputPrimaryDataset = 'eecRun3' + SAMPLES[sample][3]
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/jvelazqu/eec/' + name
    config.Data.publication = False
    config.section_('Site')
    config.Site.whitelist = ['T2_US_Vanderbilt']
    config.Site.storageSite = 'T2_CH_CERN'
    config.section_('Debug')
    config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
    return config
