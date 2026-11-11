#!/usr/bin/env python3
"""Mergeable AN7.5.1 batch reductions, reusing pair_counts physics unchanged.

Training persists raw sufficient counts, a holdout jet cache, and an exact raw
identity ledger. Normalization is performed only after global count reduction.
Closure moments reproduce the existing absolute/normalized-shape ratio errors.
"""
import argparse
import gzip
import itertools
import json
import math
import sqlite3
import struct
from collections import Counter, defaultdict
from pathlib import Path
from preparation_contract import validate_contract
import pair_counts as p
from plot_closure import FIELDS, EVENT_FIELDS, normalized_shapes

SCHEMA='run3-an751-batch-v1'
COUNT_FIELDS=tuple(p.count_cell())
YIELD_FIELDS=tuple(p.yield_cell())
IDENTITY=struct.Struct('<QQQB')


def read_json(path):
    path=Path(path)
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as stream:
        return json.load(stream)


def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.tmp')
    with (gzip.open(temporary,'wt',compresslevel=6) if path.suffix=='.gz' else temporary.open('w')) as stream:
        json.dump(data,stream,separators=(',',':'),allow_nan=False)
        stream.write('\n')
    temporary.replace(path)


def events(path):
    path=Path(path)
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as stream:
        for line in stream:
            if line.strip():yield json.loads(line)


def add(destination,source,fields=None):
    for key in source if fields is None else fields:
        destination[key]=destination.get(key,0)+source.get(key,0)


def rows(mapping):
    return [{'key':list(key),**value} for key,value in sorted(mapping.items(),key=lambda item:str(item[0]))]


def unrows(values):return {tuple(row['key']):{k:v for k,v in row.items() if k!='key'} for row in values}


def raw_counts(mapping):
    cells={};controls={}
    for key,value in mapping.items():
        cells[key]={field:value[field] for field in COUNT_FIELDS}
        control={field:value['control'][field] for field in COUNT_FIELDS}
        if key[:-1] in controls and controls[key[:-1]]!=control:
            raise ValueError('Inconsistent repeated normalization control')
        controls[key[:-1]]=control
    return {'cells':rows(cells),'controls':rows(controls)}


def train(input_path,prefix,raw_ledger=None,metadata_path=None,task=None):
    """One source file; no per-file normalized correction is persisted/applied."""
    input_path=Path(input_path);prefix=Path(prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    if Path(str(prefix)+'.receipt.json').exists():raise ValueError('Completed shard already exists')
    p.validate(p.JsonlEvents(input_path))
    table,independent,coverage=p.derive_tables(p.JsonlEvents(input_path),min_pairs=1)
    cache=Path(str(prefix)+'.holdout.jsonl.gz');identities=Path(str(prefix)+'.identities.bin.gz')
    selected={};numbers=Counter(raw=0,selected=0,train=0,holdout=0);namespaces=set();systems=set()
    with gzip.open(cache,'wt',compresslevel=6) as stream:
        for event in p.JsonlEvents(input_path):
            key=p.physical_key(event);f=p.fold(event);entry=event.get('source_entry',numbers['selected'])
            if entry in selected:raise ValueError('Repeated selected source entry')
            selected[entry]=(key,f);numbers['selected']+=1;numbers['train' if f==0 else 'holdout']+=1
            namespaces.add(key[0]);systems.add(event['sample'])
            if f==1:
                compact={k:v for k,v in event.items() if k not in ('gen','reco')}
                stream.write(json.dumps(compact,separators=(',',':'),allow_nan=False)+'\n')
    if not namespaces and task:namespaces.add(task['sample_namespace']);systems.add(task['system'])
    if len(namespaces)!=1 or len(systems)!=1:raise ValueError('Shard must contain one namespace/system')
    namespace=next(iter(namespaces));system=next(iter(systems));selected_seen=set()
    with gzip.open(identities,'wb',compresslevel=6) as stream:
        source=events(raw_ledger) if raw_ledger else ({'source_entry':entry,'id':list(key[1:])} for entry,(key,f) in sorted(selected.items()))
        for row in source:
            entry=row['source_entry'];key=(namespace,*(int(x) for x in row['id']));marker=0
            if entry in selected:
                wanted,f=selected[entry]
                if key!=wanted or entry in selected_seen:raise ValueError('Raw/selected identity mismatch')
                selected_seen.add(entry);marker=f+1
            stream.write(IDENTITY.pack(*key[1:],marker));numbers['raw']+=1
    if selected_seen!=set(selected):raise ValueError('Selected event absent from raw ledger')
    if task and numbers['raw']!=task['entry_limit']:raise ValueError('Raw coverage differs from full-file task')
    if task and (task['system']!=system or task['sample_namespace']!=namespace):raise ValueError('Task namespace/system mismatch')
    metadata=read_json(metadata_path) if metadata_path else {}
    if metadata.get('tracking_contract') is not None:
        validate_contract(metadata['tracking_contract'])
    elif task and task.get('preparation_version'):
        raise ValueError('New campaign requires corrected, versioned native export')
    result={'batch_schema':SCHEMA,'correction_version':p.VERSION,'system':system,'sample_namespace':namespace,
            'event_counts':dict(numbers),'derivation_coverage':coverage,
            'differential':raw_counts(table),'jet_independent':raw_counts(independent),
            'input_jsonl_sha256':p.file_digest(input_path),'native_metadata':metadata,'task':task,
            'raw_identity_scope':'ALL_SCANNED_EVENTS' if raw_ledger else 'SELECTED_ONLY_TEST_FIXTURE',
            'physics_core_sha256':p.file_digest(Path(p.__file__)),
            'cache_policy':'All holdout events retained, including zero-jet events; original event metadata and all jet/track fields retained; event-level auxiliary gen/reco arrays consumed only by training and removed afterward.'}
    counts=Path(str(prefix)+'.counts.json.gz');write_json(counts,result)
    receipt={'batch_schema':SCHEMA,'system':system,'sample_namespace':namespace,'event_counts':dict(numbers),'task':task,
             'files':{kind:{'name':path.name,'sha256':p.file_digest(path),'bytes':path.stat().st_size}
                      for kind,path in [('counts',counts),('holdout',cache),('identities',identities)]}}
    write_json(Path(str(prefix)+'.receipt.json'),receipt)
    return receipt


def checked_file(directory,receipt,kind):
    row=receipt['files'][kind];path=Path(directory)/row['name']
    if path.stat().st_size!=row['bytes'] or p.file_digest(path)!=row['sha256']:raise ValueError('Changed shard artifact: '+str(path))
    return path


def merge_training(manifest,directory,system,output,identity_db,min_pairs=10):
    manifest=Path(manifest);directory=Path(directory);tasks=[t for t in read_json(manifest)['tasks'] if t['system']==system]
    if not tasks:raise ValueError('No expected system tasks')
    identity_db=Path(identity_db)
    if identity_db.exists():raise ValueError('Identity database must be a new scratch file')
    connection=sqlite3.connect(str(identity_db));connection.execute('PRAGMA journal_mode=OFF');connection.execute('PRAGMA synchronous=OFF');connection.execute('PRAGMA cache_size=-32768');connection.execute('CREATE TABLE identities (namespace TEXT NOT NULL, identity BLOB NOT NULL, PRIMARY KEY(namespace,identity)) WITHOUT ROWID')
    mappings={family:{'cells':defaultdict(p.count_cell),'controls':defaultdict(p.count_cell)} for family in ['differential','jet_independent']}
    numbers=Counter(raw=0,selected=0,train=0,holdout=0);coverage=Counter();shards=[];core_hash=None;native_signature=None
    try:
        for task in tasks:
            prefix=f"{system}_{task['file_index']:03d}";receipt_path=directory/(prefix+'.receipt.json');receipt=read_json(receipt_path)
            if receipt['task']!=task:raise ValueError('Shard task does not match frozen manifest')
            count_path=checked_file(directory,receipt,'counts');counts=read_json(count_path)
            if counts['raw_identity_scope']!='ALL_SCANNED_EVENTS':raise ValueError('Full raw identity ledger required')
            if counts['batch_schema']!=SCHEMA or counts['correction_version']!=p.VERSION:raise ValueError('Incompatible shard physics version')
            core_hash=core_hash or counts['physics_core_sha256']
            if counts['physics_core_sha256']!=core_hash:raise ValueError('Mixed physics-core versions')
            signature={k:counts['native_metadata'].get(k) for k in ('analyzer_sha256','native_binary_sha256','payload_sha256','tool_sha256','event_weight_policy','jet_policy','tracking_contract')}
            if native_signature is None:native_signature=signature
            if signature!=native_signature:raise ValueError('Mixed native source/binary/calibration/policy provenance')
            ledger=checked_file(directory,receipt,'identities');cache=checked_file(directory,receipt,'holdout')
            audit=Counter(raw=0,selected=0,train=0,holdout=0)
            with gzip.open(ledger,'rb') as stream:
                batch=[]
                while True:
                    record=stream.read(IDENTITY.size)
                    if not record:break
                    if len(record)!=IDENTITY.size:raise ValueError('Truncated identity record')
                    run,lumi,event,marker=IDENTITY.unpack(record)
                    if marker not in (0,1,2):raise ValueError('Invalid identity selection marker')
                    audit['raw']+=1
                    if marker:
                        identity_event={'sample':system,'sample_namespace':task['sample_namespace'],'id':[run,lumi,event]}
                        if p.fold(identity_event)!=marker-1:raise ValueError('Identity train/holdout leakage')
                        audit['selected']+=1;audit['train' if marker==1 else 'holdout']+=1
                    batch.append((task['sample_namespace'],record[:24]))
                    if len(batch)>=10000:
                        connection.executemany('INSERT INTO identities VALUES (?,?)',batch);batch.clear()
                if batch:connection.executemany('INSERT INTO identities VALUES (?,?)',batch)
            connection.commit()
            if dict(audit)!=counts['event_counts'] or counts['event_counts']!=receipt['event_counts'] or audit['raw']!=task['entry_limit']:raise ValueError('Shard event coverage mismatch')
            numbers.update(audit);coverage.update(counts['derivation_coverage'])
            for family in mappings:
                for kind in ['cells','controls']:
                    for key,row in unrows(counts[family][kind]).items():add(mappings[family][kind][key],row,COUNT_FIELDS)
            shards.append({'prefix':prefix,'receipt':str(receipt_path.resolve()),'receipt_sha256':p.file_digest(receipt_path),
                           'holdout':str(cache.resolve()),'holdout_sha256':receipt['files']['holdout']['sha256'],
                           'event_counts':dict(audit),'input':task['input'],'input_identity':task.get('input_identity')})
        unique=connection.execute('SELECT COUNT(*) FROM identities').fetchone()[0]
        if unique!=numbers['raw']:raise ValueError('Incomplete global identity audit')
    finally:connection.close()
    result=p.run([],min_pairs=min_pairs)
    result['map']=[{'key':list(k),**v} for k,v in p.finalize_table(mappings['differential']['cells'],mappings['differential']['controls'],min_pairs).items()]
    result['jet_independent_map']=[{'key':list(k),**v} for k,v in p.finalize_table(mappings['jet_independent']['cells'],mappings['jet_independent']['controls'],min_pairs).items()]
    support={}
    for family,field in [('differential','map'),('jet_independent','jet_independent_map')]:
        for row in result[field]:
            key=tuple(row['key'][:-1]);support[(family,*key)]={'family':family,'key':list(key),**row['control'],'single_over_gen':row['large_angle_single_over_gen']}
    result['control']['support']=list(support.values());result['derivation_coverage']=dict(coverage)
    result.update(batch_schema=SCHEMA,system=system,event_counts=dict(numbers),shards=shards,
                  identity_audit={'status':'PASS','raw_unique_physical_ids':unique,'duplicate_policy':'SQLite exact namespace+uint64 run/lumi/event primary key; any duplicate aborts'},
                  input_path=str(manifest.resolve()),input_sha256=p.file_digest(manifest),physics_core_sha256=core_hash,native_signature=native_signature,
                  input_files=[{'path':t['input'],'identity':t.get('input_identity')} for t in tasks],
                  batch_status='GLOBAL_COUNTS_NORMALIZED_HOLDOUT_NOT_YET_EVALUATED')
    if native_signature and native_signature.get('tracking_contract'):
        contract = validate_contract(native_signature['tracking_contract'])
        result['tracking_contracts'] = {system: contract}
    write_json(output,result);return result


def accumulate_vectors(vectors,absolute,moments,totals,jets):
    selected_bins={i for i in range(len(p.DR)-1) if .008<=p.DR[i] and p.DR[i+1]<=.4}
    for vector in vectors:
        groups=defaultdict(dict)
        for jet in vector['jets']:
            key=tuple(jet['key']);add(jets[key],jet,('weighted_jets','jets'));jets[key]['events']+=1
        for row in vector['channels']:
            key=tuple(row['key']);values={f:row.get(EVENT_FIELDS.get(f,f),0.) for f in FIELDS};truth=values['gen']
            for field,value in values.items():
                absolute[key][field+'2']+=value*value;absolute[key][field+'_gen']+=value*truth
            if key[-1] in selected_bins:groups[key[:-1]][key[-1]]=values
        for channel,bins in groups.items():
            event_total={f:sum(row[f] for row in bins.values()) for f in FIELDS}
            for field in FIELDS:
                totals[channel][field+'2']+=event_total[field]**2;totals[channel][field+'_gen']+=event_total[field]*event_total['gen']
            for ibin,values in bins.items():
                key=(*channel,ibin)
                for field,value in values.items():
                    moments[key][field+'_self_total']+=value*event_total[field]
                    moments[key][field+'_truth_total']+=value*event_total['gen']
                    moments[key][field+'_total_truth_bin']+=event_total[field]*values['gen']


def validate_cache(cache,training,output,chunk_size=128):
    training=Path(training);result=read_json(training)
    if result.get('identity_audit',{}).get('status')!='PASS':raise ValueError('Global identity audit must pass before holdout use')
    cache=Path(cache);matching=[s for s in result['shards'] if Path(s['holdout']).name==cache.name]
    if len(matching)!=1 or p.file_digest(cache)!=matching[0]['holdout_sha256']:raise ValueError('Holdout cache not bound to global training table')
    table=unrows(result['map']);independent=unrows(result['jet_independent_map'])
    sums=defaultdict(p.yield_cell);usage={};jets=defaultdict(lambda:{'weighted_jets':0.,'jets':0,'events':0})
    absolute=defaultdict(lambda:defaultdict(float));moments=defaultdict(lambda:defaultdict(float));totals=defaultdict(lambda:defaultdict(float));number=0
    source=iter(events(cache))
    while True:
        chunk=list(itertools.islice(source,chunk_size))
        if not chunk:break
        p.validate(chunk)
        if any(p.fold(e)!=1 for e in chunk):raise ValueError('Training event leaked into holdout cache')
        detail={};spectra,vectors=p.closure(chunk,table,independent,complete_grid=False,details=detail);number+=len(chunk)
        for row in spectra:add(sums[tuple(row['key'])],row,YIELD_FIELDS)
        for row in detail['correction_application']:
            key=tuple(row['key']);target=usage.setdefault(key,{'reco_pairs':0,'single_n1':0.,'single_n2':0.,'correction_source':row['correction_source'],'correction':row['correction']})
            if target['correction_source']!=row['correction_source'] or target['correction']!=row['correction']:raise ValueError('Correction changed within fixed table')
            add(target,row,('reco_pairs','single_n1','single_n2'))
        accumulate_vectors(vectors,absolute,moments,totals,jets)
    if number!=matching[0]['event_counts']['holdout']:raise ValueError('Holdout cache coverage mismatch')
    shard={'batch_schema':SCHEMA,'prefix':matching[0]['prefix'],'training_sha256':p.file_digest(training),
           'holdout_sha256':p.file_digest(cache),'holdout_events':number,'spectra_sums':rows(sums),'jet_event_sums':rows(jets),
           'correction_application':rows(usage),'absolute_moments':rows(absolute),'shape_bin_moments':rows(moments),'shape_total_moments':rows(totals)}
    write_json(output,shard);return shard


def finalize_spectra(sums,system):
    for jet,n,selection,dr in itertools.product(range(len(p.JET)-1),(1,2),p.SELECTIONS[system],range(len(p.DR)-1)):
        sums[(system,'R04',jet,n,selection,dr)]
    spectra=[];totals=defaultdict(lambda:{'gen':0.,'raw':0.,'single':0.,'pair':0.,'complete':True})
    for key,row in sorted(sums.items()):
        populated=row['reco_pairs']>0 or row['gen_pairs']>0
        item={'key':list(key),**row,'pair':row['pair_supported'] if row['masked_pairs']==0 and populated else None,
              'pair_jet_independent':row['pair_jet_independent_supported'] if row['independent_masked_pairs']==0 and populated else None,
              'supported_fraction_single':(row['single']-row['masked_single'])/row['single'] if row['single']>0 else None,
              'status':'EMPTY' if not populated else ('UNSUPPORTED_MAP_CELLS' if row['masked_pairs'] else 'MEASURED_HOLDOUT')}
        for name in ('raw','single','pair','pair_jet_independent'):item[name+'_over_gen']=item[name]/row['gen'] if item[name] is not None and row['gen']>0 else None
        spectra.append(item);total=totals[key[:-1]]
        for name in ('gen','raw','single'):total[name]+=row[name]
        if populated:
            if item['pair'] is None:total['complete']=False
            else:total['pair']+=item['pair']
    for row in spectra:
        total=totals[tuple(row['key'][:-1])];row['shape']={f:row[f]/total[f] if total[f]>0 else None for f in ('gen','raw','single')}
        row['shape']['pair']=row['pair']/total['pair'] if row['pair'] is not None and total['complete'] and total['pair']>0 else None
    return spectra


def shape_statistics(result,absolute,moments,totals):
    normalized,_,diagnostics=normalized_shapes(result);shapes={tuple(row['key']):row for row in normalized['holdout_spectra']};integrals={tuple(row['channel']):row['integrals'] for row in diagnostics};stats={}
    roundoff=[]
    for key,row in shapes.items():
        stat={};total=integrals[key[:-1]];a=absolute.get(key,{});m=moments.get(key,{});t=totals.get(key[:-1],{})
        for field in FIELDS:
            sf=row[field] or 0.;sg=row['gen'] or 0.;tf=total[field];tg=total['gen']
            terms=(a.get(field+'2',0.),-2*sf*m.get(field+'_self_total',0.),sf*sf*t.get(field+'2',0.))
            numerator=sum(terms)
            if tf and numerator<0:
                tolerance=1e-12*max(1.,sum(abs(v) for v in terms))
                if numerator < -tolerance:raise ValueError('Negative shape variance beyond documented floating-point roundoff')
                roundoff.append({'key':list(key),'field':field,'negative_numerator':numerator,'roundoff_tolerance':tolerance})
            variance=numerator/(tf*tf) if tf else 0.
            cross=(a.get(field+'_gen',0.)-sg*m.get(field+'_truth_total',0.)-sf*m.get(field+'_total_truth_bin',0.)+sf*sg*t.get(field+'_gen',0.))/(tf*tg) if tf and tg else 0.
            stat[field+'2']=max(0.,variance);stat[field+'_gen']=cross
            # Truth self-covariance is exactly its variance. Reuse the checked
            # result so the equivalent cancellation does not round differently.
            if field=='gen':stat['gen_gen']=stat['gen2']
        stats[key]=stat
    result['shape_moment_roundoff']={'policy':'Only negative variance numerators within1e-12 of summed absolute cancellation terms are set to zero; larger negative values fail.','cells':roundoff}
    return stats


def merge_validation(training,directory,output):
    training=Path(training);result=read_json(training);digest=p.file_digest(training);directory=Path(directory)
    sums=defaultdict(p.yield_cell);usage={};jets=defaultdict(lambda:{'weighted_jets':0.,'jets':0,'events':0})
    absolute=defaultdict(lambda:defaultdict(float));moments=defaultdict(lambda:defaultdict(float));totals=defaultdict(lambda:defaultdict(float));number=0
    for expected in result['shards']:
        path=directory/(expected['prefix']+'.validation.json.gz');shard=read_json(path)
        if shard['prefix']!=expected['prefix'] or shard['training_sha256']!=digest or shard['holdout_sha256']!=expected['holdout_sha256'] or shard['holdout_events']!=expected['event_counts']['holdout']:raise ValueError('Validation shard provenance/coverage mismatch')
        number+=shard['holdout_events']
        for name,destination in [('spectra_sums',sums),('jet_event_sums',jets),('absolute_moments',absolute),('shape_bin_moments',moments),('shape_total_moments',totals)]:
            for key,row in unrows(shard[name]).items():add(destination[key],row)
        for key,row in unrows(shard['correction_application']).items():
            target=usage.setdefault(key,{'reco_pairs':0,'single_n1':0.,'single_n2':0.,'correction_source':row['correction_source'],'correction':row['correction']})
            if target['correction_source']!=row['correction_source'] or target['correction']!=row['correction']:raise ValueError('Inconsistent correction use')
            add(target,row,('reco_pairs','single_n1','single_n2'))
    if number!=result['event_counts']['holdout']:raise ValueError('Global holdout coverage mismatch')
    result['holdout_spectra']=finalize_spectra(sums,result['system']);result['jet_event_sums']=rows(jets);result['correction_application']=rows(usage)
    result['event_statistics']=rows(absolute);result['shape_statistics']=rows(shape_statistics(result,absolute,moments,totals))
    result['unsupported_holdout_bins']=[{k:row[k] for k in ('key','masked_pairs','masked_single','single','supported_fraction_single')} for row in result['holdout_spectra'] if row['masked_pairs']]
    result['covariance_status']='Exact existing absolute-bin and normalized-shape physical-event second moments; conditional on fixed training map and jet denominator. Compact holdout jet caches retain inputs for later joint covariance. Finite-training-map uncertainty not evaluated.'
    result['batch_status']='ALL_EXPECTED_SHARDS_TRAINED_AND_HELD_OUT';result['training_sha256']=digest
    write_json(output,result);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    t=sub.add_parser('train');t.add_argument('input',type=Path);t.add_argument('prefix',type=Path);t.add_argument('--raw-ledger',type=Path,required=True);t.add_argument('--metadata',type=Path,required=True);t.add_argument('--task',type=Path,required=True)
    m=sub.add_parser('merge-training');m.add_argument('manifest',type=Path);m.add_argument('directory',type=Path);m.add_argument('output',type=Path);m.add_argument('--system',choices=['OO','pp'],required=True);m.add_argument('--identity-db',type=Path,required=True);m.add_argument('--min-pairs',type=int,default=10)
    v=sub.add_parser('validate');v.add_argument('cache',type=Path);v.add_argument('training',type=Path);v.add_argument('output',type=Path)
    f=sub.add_parser('merge-validation');f.add_argument('training',type=Path);f.add_argument('directory',type=Path);f.add_argument('output',type=Path)
    a=parser.parse_args()
    if a.command=='train':r=train(a.input,a.prefix,a.raw_ledger,a.metadata,read_json(a.task))
    elif a.command=='merge-training':r=merge_training(a.manifest,a.directory,a.system,a.output,a.identity_db,a.min_pairs)
    elif a.command=='validate':r=validate_cache(a.cache,a.training,a.output)
    else:r=merge_validation(a.training,a.directory,a.output)
    print(json.dumps({'command':a.command,'event_counts':r.get('event_counts'),'holdout_events':r.get('holdout_events'),'batch_status':r.get('batch_status'),'schema':SCHEMA}))
if __name__=='__main__':main()
