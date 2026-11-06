#!/usr/bin/env python3
"""Validate, print or run a bounded native response pilot. Never unfolds."""
import argparse
import hashlib
import json
import shlex
import subprocess
import sys
from pathlib import Path
from pilot_contract import read_card, validate_card, validate_product_metadata

REPO=Path(__file__).resolve().parents[2]
DEFAULT=REPO/'validation/run3Unfolding/results/2026-11-02-lowpt-preparation/pilot_inputs.json'

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def check(manifest):
    for path,expected in manifest['source_sha256'].items():
        if digest(path)!=expected: raise ValueError('Changed source/calibration; audit and refreeze: '+path)
    if digest(manifest['full_mc_inventory_manifest'])!=manifest['full_mc_inventory_sha256']:
        raise ValueError('Changed full-MC inventory manifest')
    for task in manifest['tasks']:
        card=REPO/task['card']
        if digest(card)!=task['card_sha256']: raise ValueError('Changed pilot card: '+task['id'])
        if digest(REPO/task['file_list'])!=task['file_list_sha256']: raise ValueError('Changed input list')
        if (REPO/task['file_list']).read_text().splitlines()!=[task['input']]:
            raise ValueError('Pilot requires exactly its one frozen input file')
        stat=Path(task['input']).stat()
        identity=task['input_identity']
        if stat.st_size!=identity['size_bytes'] or stat.st_mtime_ns!=identity['mtime_ns']:
            raise ValueError('Changed raw forest: '+task['input'])
        validate_card(read_card(card),task['system'],task['power'],task['activity_edges'])
    return True

def commands(task,folder,mask):
    raw=folder/(task['id']+'_raw.root')
    projected=folder/(task['id']+'_processed.root')
    native=[str(REPO/'eecAnalysis'),str(REPO/task['file_list']),str(REPO/task['card']),str(raw),'0','true']
    projection=[sys.executable,str(REPO/'validation/run3Unfolding/project_response.py'),str(raw),str(projected),'--power',str(task['power']),'--mask',str(mask),'--scratch',str(folder/(task['id']+'_projection_build'))]
    return native,projection,raw,projected

def embedded_card(f):
    directory=f.GetDirectory('JCard')
    if not directory: raise ValueError('Missing embedded JCard')
    actual={}
    for key in directory.GetListOfKeys():
        obj=key.ReadObj()
        if obj.InheritsFrom('TVectorT<float>') or obj.ClassName().startswith('TVectorT'):
            actual[key.GetName()]=[float(obj[i]) for i in range(obj.GetLwb(),obj.GetUpb()+1)]
    return actual

def raw_report(path,task):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f=ROOT.TFile.Open(str(path))
    if not f or f.IsZombie(): raise ValueError('Unreadable native output')
    actual=embedded_card(f)
    validate_product_metadata(read_card(REPO/task['card']),actual,task['power'])
    from array import array
    histograms={}
    for name in ('jetPtUnfoldingMeasured','jetPtUnfoldingTruth','jetPtUnfoldingResponse',
                 'oneDimensionalJetPtUnfoldingMeasured','oneDimensionalJetPtUnfoldingTruth',
                 'oneDimensionalJetPtUnfoldingResponse','jetPtUnfoldingCovariance'):
        h=f.Get(name)
        if not h or not h.InheritsFrom('THnSparse'): raise ValueError('Missing response output: '+name)
        bins=array('i',[0]*h.GetNdimensions())
        totals=dict(sumw_all=0.,sumw_visible=0.,sumw_flow=0.,sumw_activity_flow=0.,activity_flow_cells=0,negative_cells=0,nonfinite_cells=0)
        activity_axis=(2 if name.endswith('Response') else 1) if name.startswith('oneDimensional') else (3 if name.endswith(('Response','Covariance')) else 2)
        import math
        for b in range(h.GetNbins()):
            weight=float(h.GetBinContent(b,bins))
            if not math.isfinite(weight): totals['nonfinite_cells']+=1; continue
            totals['negative_cells']+=int(weight<0)
            flow=any(bins[a]<1 or bins[a]>h.GetAxis(a).GetNbins() for a in range(h.GetNdimensions()))
            activity_flow=bins[activity_axis]<1 or bins[activity_axis]>h.GetAxis(activity_axis).GetNbins()
            if activity_flow and weight!=0:
                totals['sumw_activity_flow']+=weight
                totals['activity_flow_cells']+=1
            totals['sumw_all']+=weight
            totals['sumw_flow' if flow else 'sumw_visible']+=weight
        histograms[name]=dict(entries=h.GetEntries(),**totals)
        if totals['nonfinite_cells']: raise ValueError('Nonfinite response cells: '+name)
        if task['channel'] in ('MB','ref') and totals['activity_flow_cells']:
            raise ValueError('Inclusive accepted activity outside storage: expand/account bins, never impose a centrality cut')
    f.Close()
    if histograms['jetPtUnfoldingResponse']['entries']<=0:
        raise ValueError('No EEC migration entries: channel unsupported by this cap')
    return dict(status='BOUNDED_RESPONSE_ENGINEERING_ONLY',task=task['id'],histograms=histograms,
        covariance='Within-jet second moments only; full event covariance unvalidated',
        flow='Includes other activity channels for single-channel cards; inspect axes before interpreting as lost jets',
        calibration='Pair corrections disabled for pure generator response; full detector closure pending',
        normalization='Inherited generator factor1 is unverified; not an accepted cross-section normalization')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('check','commands','run'))
    parser.add_argument('--manifest',type=Path,default=DEFAULT)
    parser.add_argument('--task',action='append',default=[],help='One or more exact task IDs; otherwise firstOO and pp, both powers')
    parser.add_argument('--output-dir',type=Path,help='New output folder, outside prepared evidence recommended')
    args=parser.parse_args()
    manifest=json.loads(args.manifest.read_text())
    check(manifest)
    tasks=[t for t in manifest['tasks'] if (t['id'] in args.task if args.task else t['first_execution'])]
    if args.task and set(args.task)!={t['id'] for t in tasks}: parser.error('Unknown task ID')
    if args.action=='check':
        print('PASS:12 dedicated cards, exact serialized physical edges, powers, selections and frozen inputs. Physics readiness remains conditional.')
        return
    if args.action=='run' and args.output_dir is None: parser.error('run requires a new --output-dir')
    folder=(args.output_dir or Path('/tmp/unfolding-response-pilot')).resolve()
    if args.action=='commands':
        print('# Prepared commands only. No iteration count and no unfolding call.')
        print('cd '+shlex.quote(str(REPO)))
        for task in tasks:
            native,projection,_,_=commands(task,folder,manifest['projection_mask'])
            print(shlex.join(native))
            print(shlex.join(projection))
        return
    folder.mkdir(parents=True,exist_ok=False)
    for task in tasks:
        native,projection,raw,projected=commands(task,folder,manifest['projection_mask'])
        record=dict(task=task,commands=[native,projection],status='RUNNING')
        record_path=folder/(task['id']+'.run.json')
        record_path.write_text(json.dumps(record,indent=2)+'\n')
        try:
            with (folder/(task['id']+'.analysis.log')).open('w') as log:
                subprocess.run(native,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True)
            report=raw_report(raw,task)
            with (folder/(task['id']+'.projection.log')).open('w') as log:
                subprocess.run(projection,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=True)
            import ROOT
            f=ROOT.TFile.Open(str(projected))
            if not f or f.IsZombie() or any(not f.GetDirectory(n) for n in ('jetPtUnfoldingMeasured','jetPtUnfoldingTruth','jetPtUnfoldingResponse')):
                raise ValueError('Projection lacks required response directories')
            validate_product_metadata(read_card(REPO/task['card']),embedded_card(f),task['power'])
            f.Close()
            record.update(status='COMPLETE_ENGINEERING_ONLY',qa=report,output_sha256={str(p):digest(p) for p in (raw,projected)})
        except Exception as error:
            record.update(status='FAILED',error=str(error))
            record_path.write_text(json.dumps(record,indent=2)+'\n')
            raise
        record_path.write_text(json.dumps(record,indent=2)+'\n')
        print('Completed bounded response/projection:',task['id'])
    print('STOP: no unfolding or iteration selection. Inspect support, flow, weights and provenance before campaign expansion.')

if __name__=='__main__':
    try: main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print('BLOCKED:',error,file=sys.stderr)
        sys.exit(2)
