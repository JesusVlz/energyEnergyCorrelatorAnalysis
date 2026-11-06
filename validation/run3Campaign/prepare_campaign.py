#!/usr/bin/env python3
"""Freeze all discovered MC and write a reviewable Condor DAG; never submits."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

REPO=Path(__file__).resolve().parents[2]
TAG='R04_AN751_Response_allMC_20261103_v1'

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n')


def load_build_state(path):
    """Require an explicit existing runtime snapshot and verify its native binary."""
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError('Missing build-state file: '+str(path)+
                         '; prepare a batch-runtime snapshot and pass --build-state BUILD_STATE.json')
    build = json.loads(path.read_text())
    required = {'status', 'snapshot', 'binary_sha256', 'cmssw_base', 'scram_arch'}
    missing = required - build.keys()
    if missing:
        raise ValueError('Build state is missing fields: '+', '.join(sorted(missing)))
    if build['status'] != 'BUILT_BATCH_RUNTIME':
        raise ValueError('Build state must identify a BUILT_BATCH_RUNTIME snapshot')
    snapshot = Path(build['snapshot'])
    if not snapshot.is_absolute():
        snapshot = path.parent/snapshot
    snapshot = snapshot.resolve()
    binary = snapshot/'eecAnalysis'
    if not binary.is_file() or digest(binary) != build['binary_sha256']:
        raise ValueError('Native snapshot missing or changed; rebuild and supply its matching --build-state')
    build['snapshot'] = str(snapshot)
    return build


def freeze(old,tag):
    old = Path(old)
    if not old.is_file():
        raise ValueError('Missing input inventory: '+str(old)+'; prepare a full-MC inventory and pass --inventory MANIFEST.json')
    original=json.loads(old.read_text())
    discovered=set()
    for system,group in original['discovery'].items():
        discovered.update(str(p.resolve()) for p in Path(group['directory']).glob(group['pattern']))
    expected={str(Path(t['input']).resolve()) for t in original['tasks']}
    if discovered!=expected:
        raise ValueError('Inventory changed; refresh header inventory before submitting: '+str({'added':sorted(discovered-expected),'missing':sorted(expected-discovered)}))
    tasks=[]
    for task in original['tasks']:
        path=Path(task['input']);st=path.stat();identity=task['input_identity']
        if st.st_size!=identity['size_bytes'] or st.st_mtime_ns!=identity['mtime_ns']:
            raise ValueError('Input changed: '+str(path))
        if task['entry_limit']!=identity['entries'] or task['entry_limit']<=0:
            raise ValueError('Not a full-file task')
        tasks.append({k:task[k] for k in ['system','file_index','sample_namespace','input','entry_limit','input_identity']})
    counts=Counter(t['system'] for t in tasks)
    entries={s:sum(t['entry_limit'] for t in tasks if t['system']==s) for s in counts}
    return {'schema':'run3-allmc-condor-v1','tag':tag,'created_utc':datetime.now(timezone.utc).isoformat(),
            'tasks':tasks,'counts':dict(counts),'entries':entries,'source_inventory':str(old.resolve()),
            'source_inventory_sha256':digest(old),'input_verification':'Full current glob set plus per-file size and nanosecond mtime; header inventory full-file entry counts; workers recheck XRootD size/mtime and actual entries',
            'outputs':{s:'root://eoshome-j.cern.ch//eos/user/j/jvelazqu/HIN/eec-outputs/'+sample+'/'+tag for s,sample in [('OO','OOMC'),('pp','ppMC')]},
            'scope':{'radius':.4,'powers':[1,2],'reported_jet_pt':[40,80,120,160],'track_pt_min':2,'track_abs_eta_max':2.4,
                     'OO_inclusive':'All otherwise-selected MB events; NO centrality cut; actual source activity retained',
                     'OO_activity_channels':['inclusive','0-30','30-50','50-90','90-100_bookkeeping'],
                     'response_jet_edges':[5,10,20,40,80,120,160,200,300,500,5020],
                     'response_reco_threshold':40,'response_truth_threshold':5,
                     'training_validation':'Existing namespaced FNV1a64 physical-event identity modulo2, 0training/1validation',
                     'pair_normalization_deltaR':[.1,.4],'pair_min_pairs':10},
            'physics_pending':['Run3 inherited generator factor1 versus pp stored nonunity weights: cross-section normalization unresolved; no normalization guessed',
                'Upstream particle retention and producer jet-matching rule need production provenance',
                'Full detector calibration and pair correction qualification await independent all-MC closure',
                'OO90-100 activity-weight unity inherited/unvalidated',
                'Native response covariance is within-jet second moments; full event covariance unvalidated',
                'Support/guard stability is not established by metadata; no unfolding or iteration choice is run'],
            'authorization':'User explicitly requested preparation and submission of all-MC pair efficiency and response matrices on 2026-11-03'}

def package(snapshot,manifest,destination):
    for name in ('cardEECOOMCInclusiveR04.input','cardEECpp5v36TeVMCR04.input'):
        shutil.copy2(REPO/name,snapshot/name)
    for relative in ['validation/run3Pairs','validation/run3Unfolding','validation/run3Campaign']:
        target=snapshot/relative;target.mkdir(parents=True,exist_ok=True)
        for source in (REPO/relative).glob('*.py'):shutil.copy2(source,target/source.name)
    shutil.copytree(REPO/'validation/run3Unfolding/batch',snapshot/'validation/run3Unfolding/batch',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','results'))
    plot=snapshot/'plotting';plot.mkdir(exist_ok=True)
    for pattern in ('*.h','*.cxx','*.C'):
        for source in (REPO/'plotting').glob(pattern):shutil.copy2(source,plot/source.name)
    save(snapshot/'campaign_manifest.json',manifest)
    source_hashes={str(p.relative_to(snapshot)):digest(p) for p in sorted(snapshot.rglob('*')) if p.is_file() and p.suffix not in ('.o','.d','.pyc') and '__pycache__' not in p.parts}
    save(snapshot/'payload_files.json',source_hashes)
    with tarfile.open(destination,'w:gz') as archive:
        for path in sorted(snapshot.rglob('*')):
            if path.is_file() and path.suffix not in ('.o','.d','.pyc') and '__pycache__' not in path.parts:
                archive.add(path,arcname='analysis/'+str(path.relative_to(snapshot)),recursive=False)
    return {'archive_sha256':digest(destination),'archive_bytes':destination.stat().st_size,'native_sha256':digest(snapshot/'eecAnalysis'),'payload_files':source_hashes}

def prepare(args):
    state=args.state_dir.resolve()
    build=load_build_state(args.build_state or state/'build_state.json')
    # The historical DAG wrapper has a fixed runtime; plain pair-response queues
    # pass their build runtime explicitly instead. Reject an incompatible snapshot.
    if (build['cmssw_base'] != '/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_0_11'
            or build['scram_arch'] != 'el9_amd64_gcc12'):
        raise ValueError('DAG worker requires CMSSW_15_0_11 with el9_amd64_gcc12; use a matching build state')
    target=args.condor_root/'production/Run3MC'/args.tag
    if target.exists():raise ValueError('Use a new campaign directory: '+str(target))
    manifest=freeze(args.inventory,args.tag)
    save(state/'campaign_manifest.json',manifest)
    target.mkdir(parents=True)
    (target/'logs').mkdir()
    save(target/'campaign_manifest.json',manifest)
    shutil.copy2(REPO/'validation/run3Campaign/run_stage.sh',target/'run_stage.sh')
    payload=package(Path(build['snapshot']),manifest,target/'analysis_payload.tar.gz')
    save(target/'payload_provenance.json',{**payload,'build':build})
    sha=payload['archive_sha256']
    submit=f'''universe = vanilla
initialdir = {target}
executable = {target}/run_stage.sh
arguments = $(stage) $(index) $(system) {sha}
log = {target}/logs/jobs.$(ClusterId).log
output = {target}/logs/$(stage).$(system).$(index).$(ClusterId).out
error = {target}/logs/$(stage).$(system).$(index).$(ClusterId).err
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {target}/analysis_payload.tar.gz
transfer_output_files = ""
request_cpus = 1
request_memory = $(memory)MB
request_disk = $(disk)MB
requirements = (OpSysAndVer =?= "AlmaLinux9") && (Arch =?= "X86_64")
+JobFlavour = "$(flavour)"
+JobBatchName = "{args.tag}"
+Run3Campaign = "{args.tag}"
+Run3Stage = "$(stage)"
+Run3System = "$(system)"
+Run3TaskIndex = $(index)
queue 1
'''
    (target/'node.sub').write_text(submit)
    lines=[];nodes=[]
    def node(name,stage,index,system,memory,disk,category,flavour='workday'):
        lines.extend([f'JOB {name} node.sub',f'VARS {name} stage="{stage}" index="{index}" system="{system}" memory="{memory}" disk="{disk}" flavour="{flavour}"',f'RETRY {name} 2',f'CATEGORY {name} {category}'])
        nodes.append({'name':name,'stage':stage,'index':index,'system':system,'memory_mb':memory,'disk_mb':disk,'flavour':flavour})
    for system in ['OO','pp']:
        training=[];responses=[];validation=[]
        for index,task in enumerate(manifest['tasks']):
            if task['system']!=system:continue
            suffix=f"{system}_{task['file_index']:03d}"
            train='PT_'+suffix;response='R_'+suffix;valid='PV_'+suffix
            node(train,'pair-train',index,system,4000,8000,'NATIVE')
            node(response,'response',index,system,4000,8000,'NATIVE')
            node(valid,'pair-validate',index,system,2000,4000,'VALIDATION')
            training.append(train);responses.append(response);validation.append(valid)
        merge='PM_'+system;close='PC_'+system;rm='RM_'+system
        node(merge,'pair-merge',-1,system,4000,16000,'MERGE','tomorrow')
        node(close,'pair-close',-1,system,4000,10000,'MERGE','tomorrow')
        node(rm,'response-merge',-1,system,8000,64000,'MERGE','tomorrow')
        lines.extend(['PARENT '+' '.join(training)+' CHILD '+merge,
                      'PARENT '+merge+' CHILD '+' '.join(validation),
                      'PARENT '+' '.join(validation)+' CHILD '+close,
                      'PARENT '+' '.join(responses)+' CHILD '+rm])
    lines.extend(['MAXJOBS NATIVE 150','MAXJOBS VALIDATION 100','MAXJOBS MERGE 4','NODE_STATUS_FILE node_status.txt 30 ALWAYS-UPDATE'])
    (target/'campaign.dag').write_text('\n'.join(lines)+'\n')
    save(target/'nodes.json',nodes)
    merge_dir=args.condor_root/'merge/Run3MC'/args.tag
    merge_dir.mkdir(parents=True)
    (merge_dir/'README.txt').write_text('Dependent merge/projection/closure jobs are in the production DAG. Do not submit a second merge campaign.\nDAG: '+str(target/'campaign.dag')+'\n')
    counts=Counter(n['stage'] for n in nodes)
    result={'status':'PREPARED_NOT_SUBMITTED','directory':str(target),'dag':str(target/'campaign.dag'),
            'nodes':len(nodes),'node_counts':dict(counts),'forests':manifest['counts'],'events':manifest['entries'],
            'archive_sha256':sha,'archive_bytes':payload['archive_bytes'],'native_sha256':payload['native_sha256'],
            'manifest_sha256':digest(target/'campaign_manifest.json'),'outputs':manifest['outputs']}
    save(state/'submission.json',result);save(target/'preparation.json',result)
    (target/'README.txt').write_text('All available OO and ppRef MC, no event truncation. OO inclusive has no centrality cut.\n'
        'Submit once from this directory:\ncondor_submit_dag -maxjobs 200 -maxidle 100 campaign.dag\n'
        'Monitor: condor_q -constraint \'Run3Campaign == "'+args.tag+'"\'\n'
        'DAG node_status.txt records waiting dependencies, active/completed/failed nodes.\n'
        'Pair training -> global raw counts and exact identity audit -> independent holdout closure and ROOT/PDF outputs.\n'
        'Response train/validation matrices -> raw sum merge -> both-power inclusive/activity projections. No unfolding.\n'
        'See campaign_manifest.json physics_pending. Production output is not a physics qualification.\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default=TAG)
    parser.add_argument('--state-dir',type=Path,default=REPO/'validation/run3Campaign/results/2026-11-03-all-mc-v1',help='Preparation receipts directory; default is historical local state')
    parser.add_argument('--build-state',type=Path,help='Verified batch-runtime build JSON; default: STATE_DIR/build_state.json')
    parser.add_argument('--condor-root',type=Path,default=Path('/afs/cern.ch/user/j/jvelazqu/condor'))
    parser.add_argument('--inventory',type=Path,default=REPO/'validation/run3Pairs/results/2026-11-03-an751-all-mc-v1/expanded/sample_manifest.json')
    try:
        prepare(parser.parse_args())
    except (OSError, ValueError) as error:
        parser.error(str(error))
