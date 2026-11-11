#!/usr/bin/env python3
"""MC-only Run-3 R4 preflight/scan using native products and RooUnfoldBayes."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import re
import shlex
import sqlite3
import subprocess
import sys
import tempfile
import traceback
import numpy as np
from study_math import (REPORT_BINS, diagnostics, flattened_layout, forward_fold,
    interval_index, normalization_bins, response_operator, summarize_scan)
from study_root import (Product, REPO, assert_compatible, axis_edges, hmatrix,
    load_root, run_solver, sha256, variances, vector, write_matrix)


def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def resolve(base,value):
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (base/p).resolve()


def validate_channel(c):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',c['id']):
        raise ValueError('Unsafe channel ID')
    if c['system'] not in ('OO','pp') or c['power'] not in (1,2) or c.get('particle_pt',2)!=2:
        raise ValueError('Only R4 OO/ppRef, n=1,2, particle pT>2 GeV')
    a = c['activity']
    if len(a)!=2 or not all(math.isfinite(x) for x in a) or a[0]>=a[1]:
        raise ValueError('Explicit activity interval required')
    ordinary = [k for k in ('train_reco','train_gen','validation_reco','validation_gen') if c.get(k)]
    if ordinary and len(ordinary)!=4: raise ValueError('Both folds need RecoGen and GenGen')
    definition = c.get('truth_definition','response_marginal')
    if definition not in ('response_marginal','ordinary_GenGen'): raise ValueError('Unknown truth definition')
    if definition=='ordinary_GenGen' and len(ordinary)!=4: raise ValueError('Ordinary truth must not be replaced by a response marginal')
    if c.get('fake_scales'):
        scales = c['fake_scales']
        if set(scales)!={'A_to_B','B_to_A'} or any(not math.isfinite(float(x)) or float(x)<0 for x in scales.values()):
            raise ValueError('Both finite nonnegative fake exposures required, or neither')
    return c


def verify_split(c,base,hashes):
    """Check the supplied production index and existing worker ledgers, not new folds."""
    p = resolve(base,c['split_evidence'])
    evidence = json.loads(p.read_text())
    if evidence['product_sha256']!=hashes: raise ValueError('Evidence does not bind exactly these projected products')
    if not evidence['ledgers']: raise ValueError('Missing original worker ledgers')
    totals = [0,0]
    with tempfile.TemporaryDirectory(prefix='eec-fold-check-') as tmp:
        db = sqlite3.connect(str(Path(tmp)/'ids.sqlite'))
        db.execute('CREATE TABLE ids (key TEXT PRIMARY KEY, fold INTEGER NOT NULL)')
        try:
            for item in evidence['ledgers']:
                ledger = resolve(p.parent,item['path'])
                if sha256(ledger)!=item['sha256']: raise ValueError('Changed ledger: '+str(ledger))
                opener = gzip.open if ledger.suffix=='.gz' else open
                with opener(ledger,'rt') as f:
                    for row in csv.DictReader(f):
                        fold = int(row['fold'])
                        if fold not in (0,1) or not row['event_key']: raise ValueError('Bad identity/fold')
                        try: db.execute('INSERT INTO ids VALUES (?,?)',(row['event_key'],fold))
                        except sqlite3.IntegrityError as error: raise ValueError('Duplicated physical identity in ledgers') from error
                        totals[fold] += 1
            db.commit()
        finally: db.close()
    if not all(totals): raise ValueError('Empty event fold')
    return dict(status='SUPPLIED_INDEX_AND_LEDGER_CHECKED',events=totals,evidence_sha256=sha256(p),
        production_linkage='User-supplied index; retain original merge/worker receipts',
        embedding_source_reuse=evidence.get('embedding_source_reuse','NOT_PROVEN'))


def one_channel(root,c,base,out,action,region):
    products = []; data = {}
    keys = ('train','validation','train_reco','train_gen','validation_reco','validation_gen')
    paths = {k:resolve(base,c[k]) for k in keys if c.get(k)}
    hashes = {k:sha256(p) for k,p in paths.items()}
    if paths['train']==paths['validation'] or hashes['train']==hashes['validation']:
        raise ValueError('Train/validation products are identical')
    split = verify_split(c,base,hashes) if action=='scan' else dict(status='NOT_CHECKED_IN_PREFLIGHT')
    try:
        a = Product(root,paths['train'],c,'A'); products.append(a)
        b = Product(root,paths['validation'],c,'B'); products.append(b)
        assert_compatible(a.metadata,b.metadata,response_only=True)
        responses = [a.response(),b.response()]
        jets = a.jets
        for bounds in REPORT_BINS: interval_index(jets,bounds)
        n,dr = flattened_layout(axis_edges(responses[0][1].GetXaxis()),jets)
        for triplet in responses:
            for h in triplet[:2]:
                nr,er = flattened_layout(axis_edges(h.GetXaxis()),jets)
                if nr!=n or not np.allclose(er,dr,rtol=2e-6,atol=2e-8): raise ValueError('Fold/space axis mismatch')
        reports = []
        for s,(m,t,h) in enumerate(responses):
            counts = hmatrix(h)
            k,miss,fake = response_operator(counts,vector(m),vector(t))
            ix,iy = np.indices(counts.shape)
            if np.any(counts[(ix%n)!=(iy%n)]!=0): raise ValueError('Response violates the diagonal DeltaR construction')
            kp = np.divide(counts,counts.sum(axis=0)[None,:],out=np.zeros_like(counts),where=counts.sum(axis=0)[None,:]>0)
            reports.append(dict(K=diagnostics(k),Kprime=diagnostics(kp),
                missed_weight=miss.tolist(),fake_weight=fake.tolist(),
                zero_truth_columns=np.flatnonzero(vector(t)==0).tolist(),
                zero_efficiency_columns=np.flatnonzero((vector(t)>0)&(k.sum(axis=0)==0)).tolist()))
            data[f'response_{s}'],data[f'K_{s}'],data[f'fake_{s}'] = counts,k,fake
        data['jets'],data['dr'] = jets,dr
        record = dict(channel=c,files={k:str(v) for k,v in paths.items()},sha256=hashes,split=split,
            diagnostics=reports,cards=[a.metadata,b.metadata],root_version=root.gROOT.GetVersion(),
            requested_delta_r_region=region,reported_jet_bins=REPORT_BINS,
            status='PREFLIGHT_ONLY',truth_definition=c.get('truth_definition','response_marginal'),
            statistical_status='MEASURED_PER_JET_COVARIANCE; JOINT_CLOSURE_STATISTICS_PENDING',
            checklist_pending=['generator normalization','guard/binning stability','prior/model stress',
              'unregularized benchmark','response-MC statistics','formal bottom-line test','toy coverage',
              'tracking/pair/full-chain closure'])
        if action=='check':
            np.savez_compressed(out/'plot_data.npz',**data); dump(out/'manifest.json',record)
            return record
        references = []
        if c.get('truth_definition')=='ordinary_GenGen':
            for prefix,template in zip(('train','validation'),responses):
                reco = Product(root,paths[prefix+'_reco'],c,prefix+'_reco',True); products.append(reco)
                gen = Product(root,paths[prefix+'_gen'],c,prefix+'_gen',True); products.append(gen)
                assert_compatible(a.metadata,reco.metadata); assert_compatible(a.metadata,gen.metadata)
                rm,gt = reco.ordinary(template[0],1),gen.ordinary(template[1],3)
                references.append((rm,gt,reco.measured_cov(rm)))
        else:
            references = [(r[0],r[1],p.measured_cov(r[0])) for r,p in zip(responses,(a,b))]
        partial = out/'unfolding_scan.partial.root'
        result = root.TFile.Open(str(partial),'RECREATE')
        if not result or result.IsZombie(): raise RuntimeError('Cannot create scan output')
        try:
            for s,(training,held) in enumerate(((0,1),(1,0))):
                tm,tt,tc = responses[training]
                measured,truth,vm = references[held]
                if np.any((vector(measured)>0)&(vector(tm)<=0)): raise ValueError('Held-out measured bins lack training support')
                mask,_ = normalization_bins(dr,region)
                for bounds in REPORT_BINS:
                    j = interval_index(jets,bounds); ids = j*n+np.flatnonzero(mask)
                    if vector(tm)[ids].sum()<=0 or vector(tt)[ids].sum()<=0: raise ValueError('Unsupported reported jet bin: '+str(bounds))
                    if np.any((vector(truth)[ids]>0)&(data[f'K_{training}'][:,ids].sum(axis=0)<=0)):
                        raise ValueError('Reported truth bins have zero trained efficiency')
                # Trial iteration 1 is an API argument, not a validated choice.
                gate = root.EECUnfoldConfiguration((a,b)[held].card,0,0,c['power'],str(paths['train' if training==0 else 'validation']),1)
                response = root.RooUnfoldResponse(tm,tt,tc)
                result.cd()
                for h in (tm,tt,tc): h.Write(h.GetName())
                measured.Write(f'heldout_measured_{s}'); truth.Write(f'heldout_truth_{s}')
                data[f'measured_{s}'],data[f'truth_{s}'] = vector(measured),vector(truth)
                data[f'truth_var_{s}'],data[f'measured_cov_{s}'] = variances(truth),vm
                write_matrix(root,result,f'measured_cov_{s}',vm,axis_edges(measured.GetXaxis()),axis_edges(measured.GetXaxis()))
                for iteration in range(1,21):
                    u,vu,api = run_solver(root,response,measured,vm,iteration,f'unfolded_{s}_k{iteration}')
                    if u.GetNbinsX()!=truth.GetNbinsX(): raise ValueError('Result is not in truth space')
                    data[f'u_{s}_{iteration}'],data[f'cov_{s}_{iteration}'] = vector(u),vu
                    result.cd(); u.Write()
                    write_matrix(root,result,f'unfolded_cov_{s}_k{iteration}',vu,axis_edges(truth.GetXaxis()),axis_edges(truth.GetXaxis()))
                    exposure = c.get('fake_scales',{}).get('A_to_B' if s==0 else 'B_to_A')
                    if exposure is not None:
                        data[f'refold_{s}_{iteration}'] = forward_fold(data[f'K_{training}'],vector(u),data[f'fake_{training}'],float(exposure))
                    record['roounfold_api'] = api
                del gate
            rows,candidate = summarize_scan(data,jets,dr,region)
            record.update(status='SCAN_EXECUTED_DIAGNOSTIC_ONLY',candidate=candidate,solver_calls=40,
                refolding='FIXED_EXPOSURE_DIAGNOSTIC' if c.get('fake_scales') else 'NOT_RUN_NO_FAKE_EXPOSURE')
            root.TObjString(json.dumps(record,allow_nan=False)).Write('StudyProvenance')
        finally: result.Close()
        partial.rename(out/'unfolding_scan.root')
        np.savez_compressed(out/'plot_data.npz',**data)
        with (out/'chi2_scan.csv').open('w',newline='') as f:
            w = csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        dump(out/'candidate_iterations.json',candidate); dump(out/'manifest.json',record)
        return record
    finally:
        for p in reversed(products): p.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('check','scan'))
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True,help='New output directory')
    parser.add_argument('--drawing-library'); parser.add_argument('--roounfold-library',default='libRooUnfold')
    args = parser.parse_args(); source = json.loads(args.manifest.read_text())
    if source.get('schema_version')!=1: parser.error('Unsupported manifest version')
    channels = [validate_channel(dict(c)) for c in source['channels']]
    if not channels or len({c['id'] for c in channels})!=len(channels): parser.error('Empty or duplicate channels')
    args.output.mkdir(parents=True,exist_ok=False); (args.output/'logs').mkdir()
    (args.output/'commands.sh').write_text('#!/bin/sh\n'+shlex.join([sys.executable,*sys.argv])+'\n')
    try: root = load_root(args.drawing_library,args.roounfold_library)
    except Exception as error:
        dump(args.output/'environment.json',dict(status='BLOCKED',error=str(error))); raise
    status = []
    for c in channels:
        folder = args.output/c['id']; folder.mkdir()
        try:
            record = one_channel(root,c,args.manifest.resolve().parent,folder,args.action,source.get('delta_r_region',[.008,.39]))
            row = dict(channel=c['id'],status=record['status'],reason='')
        except Exception as error:
            (folder/'error.log').write_text(traceback.format_exc())
            row = dict(channel=c['id'],status='BLOCKED',reason=str(error))
        status.append(row); dump(args.output/'logs'/(c['id']+'.json'),row)
    with (args.output/'channel_status.csv').open('w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=['channel','status','reason']); w.writeheader(); w.writerows(status)
    sources = [Path(__file__),*Path(__file__).parent.glob('study_*.py')]
    git = subprocess.run(['git','-C',str(REPO),'rev-parse','HEAD'],text=True,capture_output=True)
    dirty = subprocess.run(['git','-C',str(REPO),'status','--porcelain'],text=True,capture_output=True)
    dump(args.output/'manifest.json',dict(input=source,input_sha256=sha256(args.manifest),
        source_sha256={str(p.relative_to(REPO)):sha256(p) for p in sources},git_head=git.stdout.strip() or 'UNAVAILABLE',
        dirty_state=dirty.stdout,channels=status,physics_status='PILOT_ONLY_NOT_FROZEN'))
    (args.output/'report.md').write_text('# Run-3 R4 study\n\n'+'\n'.join('- '+r['channel']+': '+r['status']+' '+r['reason'] for r in status)+
        '\n\nMeasured per-jet covariance is supplied. Truth Sumw2 and absent joint cross covariance make closure statistics diagnostic only. '
        'Response-MC statistics, prior/model, bottom-line, coverage and full-chain checks remain pending.\n')
    return 2 if any(r['status']=='BLOCKED' for r in status) else 0


if __name__=='__main__':
    try: sys.exit(main())
    except Exception as error:
        print('BLOCKED:',error,file=sys.stderr); sys.exit(2)
