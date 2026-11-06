#!/usr/bin/env python3
"""Prepare plain OOMC/ppMC Condor queues using the established site layout.

No submission. One forest/job produces pair-training inputs and response folds.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import sys
import tempfile

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'validation/run3Campaign'))
from prepare_campaign import freeze,package,digest,save,load_build_state

TAG='PairsResponse_R04_20261103_v1'
STATE=REPO/'validation/run3Campaign/results/2026-11-03-pairs-response-standard-v1'


def executable(path,text):
    path.write_text(text);path.chmod(0o755)


def submit_description(folder,sha,sample,system,tag,build,merge=False,kind=None):
    runtime=build['cmssw_base']+' '+build['scram_arch']+' '+sha
    if merge:
        executable_name='run_finish.sh'
        output='root://eoshome-j.cern.ch//eos/user/j/jvelazqu/HIN/eec-outputs/'+sample+'/'+tag+'/Final/'+kind+'_final.tar.gz'
        arguments=f'{kind} {system} {output} {runtime}'
        queue='queue 1';memory='8000';disk='32000';flavour='tomorrow'
        label=kind
    else:
        executable_name='run_job.sh'
        arguments='$(mode) $(input_url) $(output_url) $(stem) '+runtime+(' 10 ' if system=='OO' else ' 12 ')+'$(task_index)'
        queue='queue mode,input_url,output_url,stem,memory_mb,task_index from '+str(folder/'jobs.tsv')
        memory='$(memory_mb)';disk='8000';flavour='workday';label='PairResponse'
    return f'''# Uses the established {sample} payload/runtime/checksum-copy workflow.
universe = vanilla
initialdir = {folder}
executable = {folder}/{executable_name}
arguments = {arguments}
log = {folder}/logs/jobs.$(ClusterId).log
output = {folder}/logs/{label}.$(ClusterId).$(ProcId).out
error = {folder}/logs/{label}.$(ClusterId).$(ProcId).err
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {folder}/analysis_payload.tar.gz
transfer_output_files = ""
request_cpus = 1
request_memory = {memory}
request_disk = {disk}MB
requirements = (OpSysAndVer =?= "AlmaLinux9") && (Arch =?= "X86_64")
+JobFlavour = "{flavour}"
+JobBatchName = "{sample}_{tag}_{label}"
+Run3Campaign = "{tag}"
+Run3System = "{system}"
+Run3Stage = "{label}"
max_retries = 2
{queue}
'''


def prepare(args):
    site=args.condor_root
    state=args.state_dir.resolve()
    build=load_build_state(args.build_state or state/'build_state.json')
    snapshot=Path(build['snapshot'])
    references={sample:path for sample,path in [('OOMC',args.reference_oomc),('ppMC',args.reference_ppmc)] if path is not None}
    for sample in ('OOMC','ppMC'):
        for area in ('production','merge'):
            if (site/area/sample/args.tag).exists():raise ValueError('Use a new campaign directory')
    manifest=freeze(args.inventory,args.tag)
    manifest.update(schema='run3-standard-pair-response-v1',
                    layout='Established production/{OOMC,ppMC}/TAG, queue-from-TSV, one forest/job, direct POSC archive copy; separate merge queues',
                    fraction_guard='User-approved 3 Float_t epsilon endpoint validation; actual fractions and all Jet ID thresholds unchanged',
                    production_archive='outputs[system]/PairResponse/{system}_{file_index:03d}.tar.gz',
                    references={s:{'path':str(p),'jobs_sub_sha256':digest(p/'jobs.sub'),'worker_sha256':digest(p/'run_job.sh')} for s,p in references.items()})
    save(state/'campaign_manifest.json',manifest)
    target=snapshot/'condor';target.mkdir(exist_ok=True)
    for name in ('pair_response_worker.py','finalize_pair_response.py','frozen_input.py'):
        shutil.copy2(REPO/'condor'/name,target/name)
    with tempfile.TemporaryDirectory(prefix='eec-standard-payload-') as temporary:
        archive=Path(temporary)/'analysis_payload.tar.gz'
        payload=package(snapshot,manifest,archive)
        sha=payload['archive_sha256']
        records={}
        for system,sample in [('OO','OOMC'),('pp','ppMC')]:
            tasks=[(i,t) for i,t in enumerate(manifest['tasks']) if t['system']==system]
            production=site/'production'/sample/args.tag;production.mkdir(parents=True);(production/'logs').mkdir()
            merge=site/'merge'/sample/args.tag;merge.mkdir(parents=True);(merge/'logs').mkdir()
            output=Path('/eos/user/j/jvelazqu/HIN/eec-outputs')/sample/args.tag
            # Exactly the established prepareOOMC policy: create destinations before workers run.
            (output/'PairResponse').mkdir(parents=True,exist_ok=False)
            (output/'Final').mkdir()
            lines=[]
            for index,task in tasks:
                stem=f"{system}_{task['file_index']:03d}"
                input_url='root://eoshome-j.cern.ch/'+task['input'].replace('/eos/home-j/','/eos/user/j/')
                output_url=manifest['outputs'][system]+'/PairResponse/'+stem+'.tar.gz'
                lines.append(f'PairResponse {input_url} {output_url} {stem} 4000 {index}')
            (production/'jobs.tsv').write_text('\n'.join(lines)+'\n')
            shutil.copy2(REPO/'condor/runPairResponseJob.sh',production/'run_job.sh')
            shutil.copy2(REPO/'condor/runPairResponseFinish.sh',merge/'run_finish.sh')
            for folder in (production,merge):
                shutil.copy2(archive,folder/'analysis_payload.tar.gz')
                (folder/'analysis_payload.tar.gz.sha256').write_text(sha+'  analysis_payload.tar.gz\n')
                save(folder/'manifest.json',manifest)
                save(folder/'payload_provenance.json',{**payload,'build':build,'reference':manifest['references'].get(sample)})
            (production/'jobs.sub').write_text(submit_description(production,sha,sample,system,args.tag,build))
            executable(production/'submit.sh','#!/bin/bash\nset -euo pipefail\ncd '+shlex.quote(str(production))+'\nexec condor_submit jobs.sub\n')
            for kind in ('pairs','response'):
                (merge/('jobs-'+kind+'.sub')).write_text(submit_description(merge,sha,sample,system,args.tag,build,True,kind))
            executable(merge/'submit.sh','#!/bin/bash\nset -euo pipefail\ncase "${1:-}" in pairs|response) ;; *) echo "Usage: $0 pairs|response" >&2; exit 2;; esac\ncd '+shlex.quote(str(merge))+'\nexec condor_submit "jobs-$1.sub"\n')
            readme=f'''# {sample}: {args.tag}

Prepared using the established site queue layout. Optional reference: {references.get(sample, 'none supplied')}. No DAG.
One full forest per job, pair-training inputs and independent response folds together.
{len(tasks)} forests; {sum(t['entry_limit'] for _,t in tasks)} raw entries; both powers1,2.
OO inclusive has no centrality cut. Native selections/calibrations and source weights retained.
Response cards remain the dedicated RecoGen/40GeV-reco/5GeV-truth/covariance-enabled cards.

Production submission (once):
{site}/production/submit.sh {sample} {args.tag} --dry-run
{site}/production/submit.sh {sample} {args.tag}

After every production job succeeds, the finalizers require all expected archives and verify
all task identities/checksums before any global normalization or response merge:
{merge}/submit.sh pairs
{merge}/submit.sh response

Outputs: {output}
PairResponse/ contains one checksum-transferred archive per forest, including pair raw counts,
holdout cache and exact physical-event ledger; response train/validation ROOTs and full cards.
Final/ receives pairs_final.tar.gz (tables, closure ROOT/PDF, support summaries) and
response_final.tar.gz (raw independent matrices and both-power standard projections).
No unfolding or iteration selection is run. See manifest.json physics_pending.
'''
            (production/'README.md').write_text(readme);(merge/'README.md').write_text(readme)
            records[system]={'production':str(production),'merge':str(merge),'forests':len(tasks),'entries':sum(t['entry_limit'] for _,t in tasks),'output':str(output)}
    record={'status':'PREPARED_NOT_SUBMITTED','tag':args.tag,'prepared_utc':datetime.now(timezone.utc).isoformat(),
            'systems':records,'jobs':len(manifest['tasks']),'archive_sha256':sha,'native_sha256':build['binary_sha256'],
            'manifest_sha256':digest(state/'campaign_manifest.json'),'physics_pending':manifest['physics_pending']}
    save(state/'submission.json',record)
    print(json.dumps(record,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tag',default=TAG)
    p.add_argument('--state-dir',type=Path,default=STATE,help='Preparation receipts directory; default is historical local state')
    p.add_argument('--build-state',type=Path,help='Verified batch-runtime build JSON; default: STATE_DIR/build_state.json')
    p.add_argument('--inventory',type=Path,default=REPO/'validation/run3Pairs/results/2026-11-03-an751-all-mc-v1/expanded/sample_manifest.json',help='Frozen full-MC inventory from all_mc_campaign.py prepare')
    p.add_argument('--reference-oomc',type=Path,help='Optional existing OO campaign for provenance only')
    p.add_argument('--reference-ppmc',type=Path,help='Optional existing pp campaign for provenance only')
    p.add_argument('--condor-root',type=Path,default=Path('/afs/cern.ch/user/j/jvelazqu/condor'))
    try:
        prepare(p.parse_args())
    except (OSError, ValueError) as error:
        p.error(str(error))
