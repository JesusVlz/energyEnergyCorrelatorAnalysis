#!/usr/bin/env python3
"""Stage a frozen OO reserved-donor pool, run GenGen, and publish joint ledgers."""
import argparse
from collections import Counter
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time
from urllib.parse import urlsplit

MASK=(1<<64)-1
SUFFIXES=('.partners.csv','.pool.csv','.jets.csv','.contributions.csv')


def require(condition,message):
    if not condition: raise ValueError(message)


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda:source.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def settings(text):
    result={}
    for line in text.splitlines():
        fields=line.split('#',1)[0].split()
        if not fields: continue
        require(fields[0] not in result,'Duplicate card setting: '+fields[0])
        result[fields[0]]=fields[1:]
    return result


def scalar(card,key):
    require(key in card and len(card[key])==1,'Missing scalar '+key)
    value=float(card[key][0]);require(math.isfinite(value),'Nonfinite '+key)
    return value


def rewrite_card(text,pool_list,prefix):
    replacements={'MixingFileList':str(pool_list),'OOClosureOutputPrefix':str(prefix)}
    require(all(key in settings(text) for key in replacements),'Frozen closure card lacks local path settings')
    return '\n'.join((line.split()[0]+' '+replacements[line.split()[0]])
                     if line.strip() and not line.lstrip().startswith('#') and line.split()[0] in replacements else line
                     for line in text.splitlines())+'\n'


def load_plan(analysis,card_path):
    frozen=card_path.read_text();card=settings(frozen)
    for key,value in {'OOClosurePilot':1,'OOClosurePartitionPolicy':1,'DataType':10,'McCorrelationType':3,
                      'BackgroundMethods':6,'EventMatchForMixing':2,'BinnedMixing':0,'MegaSkimMixing':0,
                      'SubeventCut':2,'JetAxis':1,'JetRadius':.4,'HFEnergyShift':0,'SmearDeltaR':0,
                      'SmearEnergyWeight':0,'DisableTrackPairEfficiencyCorrection':1}.items():
        require(scalar(card,key)==value,'Reserved OO GenGen worker requires '+key+' '+str(value))
    require(card.get('WeightExponent')==['1','2'],'Worker requires both powers 1 and 2')
    manifest_path=analysis/'mixing_pool/manifest.json'
    manifest=json.loads(manifest_path.read_text())
    require(manifest.get('schema_version')==1,'Unsupported pool manifest version')
    require(manifest.get('partition_policy')==1,'Campaign must reserve the same fixed donor pool')
    require(manifest.get('holdout_policy')=='fixed_reserved_donors_exclude_overlapping_signals','Unknown global holdout policy')
    require(scalar(card,'OOClosureMaxPoolEvents')==manifest['global_event_cap'],'Pool cap differs from frozen manifest')
    require(scalar(card,'OOClosureSeed')==manifest['seed'],'Seed differs from frozen pool')
    require(scalar(card,'OOClosureMaxSignalEvents')>=manifest['signal_max_entries'],'Signal cap truncates a full source forest')
    require(manifest['global_event_cap']>0 and manifest['entries'],'Empty frozen donor prefix')
    seen=set();total=0
    for entry in manifest['entries']:
        name=entry['local_name'];require(re.fullmatch(r'pool_[0-9]+\.root',name) is not None,'Unsafe donor local filename')
        require(name not in seen,'Duplicate staged donor name');seen.add(name)
        require(entry['input_url'].startswith('root://'),'Donor URL must use XRootD')
        require(0<entry['use_entries']<=entry['source_entries'],'Invalid bounded donor prefix')
        require(entry['size_bytes']>0 and entry['root_uuid'],'Missing frozen donor identity/size')
        total+=entry['use_entries']
    require(total==manifest['global_event_cap'],'Donor prefix does not equal global cap')
    require(all(entry['use_entries']==entry['source_entries'] for entry in manifest['entries'][:-1]),'Only the last staged donor file may be partially used')
    frozen_list=analysis/manifest['frozen_pool_list']
    require(card['MixingFileList']==[manifest['frozen_pool_list']],'Frozen card pool path differs from manifest')
    require(sha256(frozen_list)==manifest['pool_list_sha256'],'Frozen pool list hash mismatch')
    listed=[line.strip() for line in frozen_list.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
    require(listed==[entry['original_path'] for entry in manifest['entries']],'Frozen donor paths/order differ from staged manifest')
    for name,field in [('source_audit.manifest.json','source_audit_sha256'),('signal_headers.manifest.json','signal_headers_sha256')]:
        require(sha256(analysis/'mixing_pool'/name)==manifest[field],'Frozen '+name+' hash mismatch')
    return frozen,card,manifest


def root_info(path,output=False):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    file=ROOT.TFile.Open(str(path))
    require(bool(file) and not file.IsZombie() and not file.TestBit(ROOT.TFile.kRecovered),'Unreadable/zombie/recovered ROOT '+str(path))
    try:
        if output:
            events=file.Get('nEvents');jets=file.Get('inclusiveJet');eec=file.Get('energyEnergyCorrelator')
            require(bool(events) and bool(jets) and bool(eec),'Missing standard output histograms')
            return {'all_events':float(events.GetBinContent(1)),'inclusive_jet_entries':float(jets.GetEntries()),'eec_entries':float(eec.GetEntries())}
        names=('hiEvtAnalyzer/HiTree','HiGenParticleAna/hi','skimanalysis/HltTree','akCs4PFJetAnalyzer/t')
        trees=[file.Get(name) for name in names];require(all(bool(tree) for tree in trees),'Missing native signal/donor forest trees')
        counts=[int(tree.GetEntries()) for tree in trees];require(len(set(counts))==1 and counts[0]>0,'Unaligned/empty source forest')
        return {'source_entries':counts[0],'root_uuid':str(file.GetUUID().AsString()),'size_bytes':Path(path).stat().st_size}
    finally:file.Close()


def check_staged(path,expected,inspect=root_info):
    require(Path(path).is_file() and Path(path).stat().st_size==expected['size_bytes'],'Staged ROOT size differs from frozen source')
    actual=inspect(path)
    for key in ('source_entries','root_uuid','size_bytes'):require(actual[key]==expected[key],'Staged ROOT '+key+' differs from frozen source')
    return actual


def csvrows(path):
    with Path(path).open(newline='') as source:
        reader=csv.DictReader(source);require(reader.fieldnames is not None,'Missing ledger header '+str(path))
        for row in reader:
            require(None not in row and None not in row.values(),'Incomplete ledger row '+str(path));yield row


def uid(value):
    require(re.fullmatch(r'[0-9]+',str(value)) is not None,'Identity must be an exact decimal uint64')
    number=int(value);require(0<=number<=MASK,'Identity outside uint64');return number


def physical(row,prefix=''):return tuple(uid(row[prefix+key]) for key in ('run','lumi','event'))


def number(value):
    result=float(value);require(math.isfinite(result),'Nonfinite ledger value');return result


def validate_ledgers(prefix,source_entries,pool_cap,seed):
    pool={};reserved_ids=set();reserved_bg=set();groups=Counter()
    for row in csvrows(str(prefix)+'.pool.csv'):
        key=uid(row['donor_key']);require(key not in pool,'Duplicate donor physical key')
        require(row['pool'] in ('A','B') and row['eligible'] in ('0','1'),'Invalid pool eligibility/group')
        bg=uid(row['background_key']);fold_key=bg or key
        require(row['pool']==('B' if (fold_key^seed)%2 else 'A'),'Pool fold violates frozen identity partition')
        require(number(row['generator_weight'])==1.,'Nonunity donor weights need a separately approved sampling study')
        pool[key]=row;reserved_ids.add(physical(row));groups[row['pool']]+=1
        if bg:reserved_bg.add(bg)
    require(len(pool)==pool_cap,'Native donor scan did not equal frozen global pool cap')
    assignments={};uses=Counter();failures=Counter();reserved_failures=0
    for row in csvrows(str(prefix)+'.partners.csv'):
        key=uid(row['signal_key']);require(key not in assignments,'Repeated physical signal assignment')
        require(row['status'] in ('matched','failed'),'Unknown event assignment status')
        bg=uid(row['signal_background_key']);reserved=physical(row) in reserved_ids or (bg and bg in reserved_bg)
        number(row['event_weight']);require(int(row['inspected']) in (0,pool_cap),'Partner search exceeds/does not follow bounded pool')
        if reserved:
            require(row['status']=='failed' and row['reason'] in ('signal_in_reserved_donor_physical_id','signal_in_reserved_donor_background'),'Reserved donor leaked into signal exposure')
            reserved_failures+=1
        if row['status']=='matched':
            require(not reserved and row['reason']=='matched','Invalid matched event reason')
            donor_keys=[]
            for partner,group in ((1,'A'),(2,'B')):
                base='donor'+str(partner)+'_';donor_key=uid(row[base+'key']);require(donor_key in pool,'Partner outside frozen donor pool')
                donor=pool[donor_key];require(donor['eligible']=='1' and donor['pool']==group,'Ineligible/wrong-pool partner')
                require(physical(row,base)==physical(donor),'Partner physical identity mismatch')
                require(abs(number(row['vz'])-number(donor['vz']))<=.5,'Partner violates fixed dz')
                require(abs(int(row['nh'])-int(donor['nh']))<=max(1.,.01*int(row['nh'])),'Partner violates fixed NH')
                donor_keys.append(donor_key);uses[donor_key]+=1
            require(donor_keys[0]!=donor_keys[1],'M1 equals M2')
        else:
            require(row['reason'] in ('no_pool_A','no_pool_B','shared_partner_background','signal_in_reserved_donor_physical_id','signal_in_reserved_donor_background'),'Unknown matching failure')
            require(all(row['donor'+str(p)+'_'+field]=='' for p in (1,2) for field in ('key','run','lumi','event','entry','vz','nh','weight')),'Partial estimator on failed event')
            failures[row['reason']]+=1
        assignments[key]=row
    require(len(assignments)<=source_entries,'More selected physical signals than source events')
    for key,row in pool.items():require(int(row['uses'])==uses[key],'Donor reuse does not equal event assignments')
    jets={};matched_jets=0
    for row in csvrows(str(prefix)+'.jets.csv'):
        key=uid(row['signal_key']);index=(key,int(row['jet_index']))
        require(key in assignments and row['status']==assignments[key]['status'],'Jet did not reuse its event assignment')
        require(index not in jets,'Repeated jet exposure');jets[index]=row
        if row['status']=='matched':matched_jets+=1
    contributions=0
    for row in csvrows(str(prefix)+'.contributions.csv'):
        index=(uid(row['signal_key']),int(row['jet_index']))
        require(index in jets and jets[index]['status']=='matched','Contribution outside common matched jet exposure')
        require(number(row['power']) in (1.,2.),'Unexpected energy power');number(row['value']);contributions+=1
    require(matched_jets or not contributions,'Nonzero EEC contributions with zero matched exposure')
    return {'source_entries':source_entries,'selected_signal_events':len(assignments),'matched_signal_events':sum(row['status']=='matched' for row in assignments.values()),
            'matched_jets':matched_jets,'all_selected_jet_rows':len(jets),'contribution_rows':contributions,'failures':dict(failures),
            'reserved_signal_failures':reserved_failures,'pool_entries':len(pool),'pool_groups':dict(groups),
            'unique_assigned_donors':len(uses),'partner_assignments':sum(uses.values()),'maximum_partner_reuse':max(uses.values(),default=0)}


def run_job(args,execute=None,inspect=root_info):
    analysis=Path(args.analysis_dir).resolve();work=Path(args.work_dir).resolve();work.mkdir(parents=True,exist_ok=True)
    card_path=analysis/args.card;frozen,card,manifest=load_plan(analysis,card_path)
    require(re.fullmatch(r'[A-Za-z0-9_.-]+',args.stem) is not None and args.stem not in ('.','..'),'Unsafe output stem')
    require(args.mode=='GenGen' and args.data_type==10,'Mixed worker supports OO GenGen only')
    require(args.input_url.startswith('root://') and args.output_url.startswith('root://') and args.output_url.endswith('.root'),'Invalid XRootD input/output URL')
    require(re.fullmatch(r'[a-fA-F0-9]{64}',args.archive_sha256) is not None,'Invalid archive SHA256')
    metadata={'status':'PREFLIGHT_PASS' if args.preflight_only else 'RUNNING','mode':'GenGen','data_type':10,'input_url':args.input_url,'output_url':args.output_url,
        'archive_sha256':args.archive_sha256,'frozen_card_sha256':sha256(card_path),'pool_manifest_sha256':sha256(analysis/'mixing_pool/manifest.json'),
        'source_audit_sha256':manifest['source_audit_sha256'],'signal_headers_sha256':manifest['signal_headers_sha256'],
        'executable_sha256':sha256(analysis/'eecAnalysis'),'partition_policy':1,'holdout_policy':manifest['holdout_policy'],
        'reserved_pool_scope':manifest.get('reserved_pool_scope'),'matching':{'nh':'stable charged HIJING pT>0.7 |eta|<2.4','dz_max_cm':.5,'dnh':'max(1,0.01 NH_signal)'},
        'seed':manifest['seed'],'pool_cap':manifest['global_event_cap'],'signal_cap':int(scalar(card,'OOClosureMaxSignalEvents')),
        'corrections':'Native GenGen particle/jet kinematics; reconstruction single-track/pair corrections unity; Run-3 reco jet veto kept.',
        'limitations':['Stored kinematic fingerprints do not certify complete production lineage.','Fixed reserved donor pool and partner reuse induce correlations; sidecars are joint-statistics input, not final covariance.'],
        'started_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'commands':[]}
    provenance_path=analysis/'campaign_provenance.json'
    require(provenance_path.is_file(),'Missing frozen campaign provenance record')
    metadata['campaign_provenance']=json.loads(provenance_path.read_text())
    metadata['campaign_provenance_sha256']=sha256(provenance_path)
    recorded=metadata['campaign_provenance'].get('executable_sha256')
    require(recorded==metadata['executable_sha256'],'Worker executable differs from frozen campaign provenance')
    if args.preflight_only:
        print(json.dumps(metadata,sort_keys=True));return metadata
    prefix=work/args.stem;output=Path(str(prefix)+'.root');run_record=Path(str(prefix)+'.run.json')
    def write_record():run_record.write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    def command(command):
        started=time.monotonic()
        if execute is None:
            result=subprocess.run(command,cwd=analysis,check=False);code=result.returncode
        else:code=execute(command,analysis)
        metadata['commands'].append({'argv':[str(item) for item in command],'exit_code':code,'seconds':time.monotonic()-started})
        require(code==0,'Worker command failed: '+str(command))
    try:
        staged=work/'mixing_pool';staged.mkdir()
        local_pool=[]
        for entry in manifest['entries']:
            local=staged/entry['local_name'];command(['xrdcp','--retry','3','--cksum','adler32',entry['input_url'],str(local)])
            check_staged(local,entry,inspect);local_pool.append(str(local))
        pool_list=work/'local_pool.txt';pool_list.write_text('\n'.join(local_pool)+'\n')
        local_signal=work/'input.root';command(['xrdcp','--retry','3','--cksum','adler32',args.input_url,str(local_signal)])
        headers=json.loads((analysis/'mixing_pool/signal_headers.manifest.json').read_text())
        filename=Path(urlsplit(args.input_url).path).name
        candidates=[row for row in headers['files'] if Path(row['path']).name==filename]
        require(len(candidates)==1,'Signal URL not uniquely present in frozen headers')
        source=check_staged(local_signal,candidates[0],inspect)
        require(source['source_entries']<=scalar(card,'OOClosureMaxSignalEvents'),'Signal global cap truncates this source forest')
        require('MaxEventsPerFile' not in card or scalar(card,'MaxEventsPerFile')<=0 or scalar(card,'MaxEventsPerFile')>=source['source_entries'],'MaxEventsPerFile truncates this full-file campaign job')
        runtime_card=work/'runtime.input';runtime_card.write_text(rewrite_card(frozen,pool_list,prefix))
        before,after=settings(frozen),settings(runtime_card.read_text())
        require({key for key in before if before[key]!=after[key]}<={'MixingFileList','OOClosureOutputPrefix'},'Worker changed frozen physics settings')
        metadata.update(source=source,source_identity={'original_path':candidates[0]['path'],'root_uuid':source['root_uuid'],'full_event_entries':source['source_entries']},
                        runtime_card_sha256=sha256(runtime_card),local_pool_list_sha256=sha256(pool_list),staged_pool=manifest['entries'])
        command([str(analysis/'eecAnalysis'),str(local_signal),str(runtime_card),str(output),'1','true'])
        ledger=validate_ledgers(prefix,source['source_entries'],manifest['global_event_cap'],manifest['seed'])
        actual=inspect(output,output=True)
        require(actual['all_events']==source['source_entries'],'ROOT all-event count does not cover the full signal forest')
        require(actual['inclusive_jet_entries']==ledger['matched_jets'],'ROOT jet exposure differs from joint ledger')
        empty=ledger['matched_jets']==0
        macro='condor/checkOOMCOutput.C('+json.dumps(str(output))+',"GenGen",10,'+('true' if empty else 'false')+')'
        command(['root','-l','-b','-q',macro])
        artifacts=[output,*[Path(str(prefix)+suffix) for suffix in SUFFIXES]]
        require(all(path.is_file() for path in artifacts),'Missing output ROOT/joint ledger sidecar')
        metadata.update(status='PASS_FULL_FILE_SOFTWARE_AND_EXPOSURE_CHECKS',ledger_validation=ledger,root_validation=actual,
                        empty_common_exposure=empty,artifacts={path.name:{'sha256':sha256(path),'size_bytes':path.stat().st_size} for path in artifacts})
        publication=[]
        for path in artifacts:
            destination=args.output_url if path==output else args.output_url[:-5]+path.name[len(args.stem):]
            command(['xrdcp','--retry','3','--cksum','adler32','--posc','--force',str(path),destination]);publication.append(destination)
        metadata.update(published_artifacts=publication,run_metadata_url=args.output_url[:-5]+'.run.json',finished_utc=dt.datetime.now(dt.timezone.utc).isoformat())
        write_record();command(['xrdcp','--retry','3','--cksum','adler32','--posc','--force',str(run_record),metadata['run_metadata_url']])
        print(json.dumps({'status':metadata['status'],'output':args.output_url,'matched_jets':ledger['matched_jets']}));return metadata
    except Exception as error:
        metadata.update(status='FAILED',error=str(error),finished_utc=dt.datetime.now(dt.timezone.utc).isoformat());write_record();raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('analysis-dir','work-dir','card','input-url','output-url','stem','archive-sha256'):parser.add_argument('--'+name,required=True)
    parser.add_argument('--mode',default='GenGen');parser.add_argument('--data-type',type=int,default=10)
    parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args();run_job(args)


if __name__=='__main__':main()
