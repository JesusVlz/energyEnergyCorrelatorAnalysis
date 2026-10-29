#!/usr/bin/env python3
"""Bounded, paired R4 calibration diagnostics; no candidate becomes nominal.

Use the legacy WTA lookup convention for both payload sets. The optional standard
coordinate variation changes only lookup coordinates, with the installed payload.
No L1 is prepended without an established producer/rho prescription. Existing
trackMax/rawpt selection is fixed; no veto or Jet ID change is made here.
"""
import argparse
import hashlib
import json
import math
import re
import shlex
import subprocess
import sys
from array import array
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ROOT

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'validation/run3Pairs'))
from export_forest import digest, select_branches

PT_EDGES=np.array([0,10,20,40,60,80,120,160,250,500.])
GROUPS={'integrated':(-1,200),'0_30':(0,60),'30_50':(60,100),'50_90':(100,180),'90_100':(180,200)}


def output_bin(pt):
    return int(np.searchsorted([40.,80.,120.,160.],pt,side='right'))


def native_match(jets,index):
    # Byte-for-byte convention of HighForestReader::GetMatchingIndex: no rematch,
    # no new wrapped-phi criterion or calibration-dependent pT-ratio veto.
    if float(jets.refpt[index])<0:return None
    for gen in range(int(jets.ngen)):
        if (abs(float(jets.refeta[index])-float(jets.geneta[gen]))<.015 and
            abs(float(jets.refphi[index])-float(jets.genphi[gen]))<.015 and
            abs(float(jets.refpt[index])-float(jets.genpt[gen]))<.03*float(jets.genpt[gen])):
            return gen
    return None


def weighted_profile(x,y,weights,edges):
    rows=[]
    for low,high in zip(edges[:-1],edges[1:]):
        keep=(x>=low)&(x<high)&np.isfinite(y)&np.isfinite(weights)
        a,b=y[keep],weights[keep];total=float(b.sum());total2=float(np.square(b).sum())
        if not len(a) or total<=0:
            rows.append(dict(low=float(low),high=float(high),count=int(len(a)),sumw=total,mean=None,width=None,neff=None));continue
        mean=float(np.dot(a,b)/total);variance=float(np.dot(np.square(a-mean),b)/total)
        rows.append(dict(low=float(low),high=float(high),count=int(len(a)),sumw=total,
                         mean=mean,width=math.sqrt(variance) if variance>=0 else None,
                         neff=total*total/total2 if total2>0 else None))
    return rows


def read_weight_coefficients():
    source=(REPO/'src/EECAnalyzer.cxx').read_text()
    oo=source[source.index('} else if(fIsOOData){'):source.index('else if(fIsPp2024Data)')]
    pp=source[source.index('else if(fIsPp2024Data)'):source.index('} else if(fIsPPbData){')]
    def params(block,name):
        text=re.search(re.escape(name)+r'->SetParameters\((.*?)\);',block,re.S).group(1)
        return [float(x.strip()) for x in text.split(',')]
    assert 'fPtHatWeight = run3MC ? 1.0 : fJetReader->GetEventWeight();' in source
    assert 'if(hiBin < 0 || hiBin >= 180) return 1.0;' in source
    return {'OO_vz':params(oo,'fVzWeightFunction'),'OO_central':params(oo,'fCentralityWeightFunctionCentral'),
            'OO_peripheral':params(oo,'fCentralityWeightFunctionPeripheral'),'pp_vz':params(pp,'fVzWeightFunction')}


def event_weight(system,vz,hi,coefficients):
    def evaluate(values,x):return sum(value*x**i for i,value in enumerate(values))
    weight=evaluate(coefficients[system+'_vz'],vz)
    # Preserve the analyzer's existing GetHiBin() negative-coordinate convention.
    hi_for_lookup=max(hi,1) if hi<0 else hi
    if system=='OO' and hi_for_lookup<180:
        weight*=evaluate(coefficients['OO_central' if hi_for_lookup<60 else 'OO_peripheral'],hi_for_lookup/2.)
    return weight


class Corrector:
    def __init__(self,files):
        self.paths=[Path(p) for p in files]
        self.levels=[ROOT.SingleJetCorrector(str(p)) for p in self.paths]
        self.bounds=[]
        for path in self.paths:
            lines=path.read_text().splitlines();header=lines[0].strip('{}').split()
            bins=int(header[0]);deps=int(header[1+bins]);names=header[2+bins:2+bins+deps]
            if bins!=1 or header[1]!='JetEta' or names!=['JetPt']:
                raise ValueError('This diagnostic requires eta/pT-only levels; no invented rho or area: '+str(path))
            self.bounds.append([(list(map(float,line.split()))[0:2],list(map(float,line.split()))[3:5])
                                for line in lines[1:] if line.strip() and not line.lstrip().startswith('#')])
    def evaluate(self,pt,eta,phi):
        # Use the same scalar level evaluator sequentially. The wrapper never reads
        # or supplies uninitialized/default rho/area for eta/pT-only payloads.
        result=pt
        for level in self.levels:
            level.SetJetPT(result);level.SetJetEta(eta);level.SetJetPhi(phi)
            result=float(level.GetCorrectedPT())
            if not math.isfinite(result) or result<=0:raise ValueError('Invalid correction output')
        clamp=False
        for eta_range,pt_range in self.bounds[0]:
            if eta_range[0]<=eta<=eta_range[1]:
                clamp=not pt_range[0]<=pt<=pt_range[1];break
        return result,clamp


def setup_library(out):
    flags=shlex.split(subprocess.check_output(['root-config','--cflags','--libs'],text=True))
    library=out/'JetCorrectorAudit.so'
    command=['g++','-shared','-fPIC','-O2',str(REPO/'src/JetCorrector.cxx'),*flags,'-o',str(library)]
    subprocess.run(command,check=True)
    assert ROOT.gSystem.Load(str(library))>=0
    ROOT.gInterpreter.Declare('#include "'+str(REPO/'src/JetCorrector.h')+'"')
    ROOT.TF1.DefaultAddToGlobalList(False)
    return command


def scan(system,path,cap,out,upstream,coefficients,data=False):
    old_l2=REPO/'jetEnergyCorrections'/('L2Relative_AK4PF_OO_V1.txt' if system=='OO' else 'Prompt24HIpp_V1_MC_L2Relative_AK4PF.txt')
    new_l2=upstream/('txt/2025OO/Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt' if system=='OO' else 'txt/2024ppRef/Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt')
    old_files=[old_l2];new_files=[new_l2]
    if data and system=='pp':
        old_files.append(REPO/'jetEnergyCorrections/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt')
        new_files.append(upstream/'txt/2024ppRef/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt')
    correctors={'installed':Corrector(old_files),'candidate_raw_input':Corrector(new_files),
                'coordinate_only':Corrector(old_files)}
    if data and system=='pp':
        correctors['residual_only_candidate']=Corrector([old_l2,new_files[-1]])
        correctors['l2_only_candidate']=Corrector([new_l2,old_files[-1]])
    f=ROOT.TFile.Open(str(path));assert f and not f.IsZombie()
    tree_name=('akCs4' if system=='OO' else 'ak4')+'PFJetAnalyzer/t'
    event,jets,skim=f.Get('hiEvtAnalyzer/HiTree'),f.Get(tree_name),f.Get('skimanalysis/HltTree')
    assert event and jets and skim
    event_branches=['run','lumi','evt','vz','hiBin']+([] if data else ['weight'])
    select_branches(event,event_branches)
    fields=['nref','rawpt','jtpt','jteta','jtphi','WTAeta','WTAphi','trackMax']
    if not data:fields+=['refpt','refeta','refphi','ngen','genpt','geneta','genphi','WTAgeneta','WTAgenphi']
    select_branches(jets,fields)
    filters=['pprimaryVertexFilter','OOphfCoincFilterPF2Th4','pileupVertexFilter'] if system=='OO' else ['pprimaryVertexFilter','ppNoScrapingFilter']
    select_branches(skim,filters)
    assert event.GetEntries()==jets.GetEntries()==skim.GetEntries()
    count=Counter();rows=[];identities=set();limit=min(cap,int(event.GetEntries()))
    for entry in range(limit):
        event.GetEntry(entry);skim.GetEntry(entry);count['read_events']+=1
        identity=(int(event.run),int(event.lumi),int(event.evt))
        if identity in identities:raise ValueError('Duplicate physical event ID')
        identities.add(identity)
        if abs(float(event.vz))>15 or not all(bool(getattr(skim,x)) for x in filters):continue
        count['event_quality_pass']+=1;jets.GetEntry(entry)
        hibin=int(event.hiBin);weight=1. if data else event_weight(system,float(event.vz),hibin,coefficients)
        if system=='OO' and not 0<=hibin<200:count['invalid_stored_hiBin_selected']+=1
        for j in range(int(jets.nref)):
            raw=float(jets.rawpt[j]);eta=float(jets.WTAeta[j]);phi=float(jets.WTAphi[j])
            if raw<=0 or abs(eta)>=1.6 or not .01<float(jets.trackMax[j])/raw<.98:continue
            count['fixed_quality_jets_before_analysis_pt']+=1
            match=None if data else native_match(jets,j)
            if not data and match is None:count['unmatched_jets']+=1;continue
            row={'id':list(identity),'source_entry':entry,'jet_index':j,'gen_index':match,'system':system,
                 'data':data,'hiBin':hibin,'centrality':hibin/2.,'event_weight':weight,
                 'stored_generator_weight':None if data else float(event.weight),
                 'rawpt':raw,'stored_jtpt':float(jets.jtpt[j]),'WTAeta':eta,'WTAphi':phi,
                 'jteta':float(jets.jteta[j]),'jtphi':float(jets.jtphi[j]),
                 'genpt':None if data else float(jets.genpt[match]),
                 'geneta':None if data else float(jets.geneta[match])}
            for name,corrector in correctors.items():
                coordinate=(row['jteta'],row['jtphi']) if name=='coordinate_only' else (eta,phi)
                value,clamped=corrector.evaluate(raw,*coordinate);row[name]=value
                count[name+'_first_level_clamped']+=int(clamped)
            row['candidate_over_installed']=row['candidate_raw_input']/row['installed']
            row['old_bin']=output_bin(row['installed']);row['candidate_bin']=output_bin(row['candidate_raw_input'])
            rows.append(row)
    uuid=str(f.GetUUID().AsString());f.Close()
    label=system+('_data' if data else '_MC')
    ledger=out/(label+'_paired_jets.jsonl')
    ledger.write_text(''.join(json.dumps(r,allow_nan=False)+'\n' for r in rows))
    return rows,{'label':label,'system':system,'data':data,'input':str(path.resolve()),'input_sha256':digest(path),
                 'input_uuid':uuid,'entry_range':[0,limit],'native_tree':tree_name,'counts':dict(count),
                 'common_sample':'fixed event quality, abs(WTAeta)<1.6, fixed trackMax/rawpt cut, native legacy ref/gen matching; BEFORE analysis jet-pT cuts; no pT-ratio rematch or new veto/ID',
                 'weights':'unity data diagnostic; existing Run-3 analyzer vz*OO-centrality, generator/cross-section factor unity for MC; stored generator branch retained separately',
                 'data_limitations':'event-quality-only data diagnostic, trigger applicability NOT_VALIDATED; no MC truth response' if data else '',
                 'payloads':{name:[{'path':str(p),'sha256':digest(p)} for p in c.paths] for name,c in correctors.items()},
                 'candidate_status':'NOT_APPROVED: upstream payload numerical comparison, raw-input hypothesis, no standalone L1; V2 includes-L1 and CHS/reconstruction applicability unresolved',
                 'coordinate_diagnostic':'same installed levels and WTA cones/sample; vary only lookup eta/phi, not a nominal change',
                 'energy_powers':{'1':'same selected jets/matches/event weights','2':'same selected jets/matches/event weights; jet calibration is power-independent; pair-weighted EEC closure NOT_RUN'},
                 'veto':{'applied':False,'event_loss_fraction':None,'before_after_maps':'before-only diagnostic; after NOT_RUN until era/map-key/ID/pT applicability is verified'},
                 'ledger':str(ledger),'ledger_sha256':digest(ledger)}


def save_fig(fig,path):
    fig.tight_layout();fig.savefig(str(path)+'.png',dpi=145);fig.savefig(str(path)+'.pdf');plt.close(fig)


def make_plots(rows,metadata,out):
    label=metadata['label'];data=metadata['data'];system=metadata['system']
    rootfile=ROOT.TFile(str(out/(label+'_diagnostics.root')),'RECREATE')
    profiles={};artifacts=[]
    if data and system=='pp':
        fig,axes=plt.subplots(2,2,figsize=(11,7.5))
        fig.suptitle('pp data: residual-only vs L2-only candidate changes; fixed quality-only events; NOT approved',fontsize=10)
        weight=np.array([r['event_weight'] for r in rows])
        variants=['installed','residual_only_candidate','l2_only_candidate','candidate_raw_input']
        for variant in variants:
            values=np.array([r[variant] for r in rows]);hist,edges=np.histogram(values,bins=np.linspace(0,250,51),weights=weight)
            axes[0,0].step(edges[:-1],hist,where='post',label=variant)
            if variant=='installed':continue
            ratios=np.array([r[variant]/r['installed'] for r in rows])
            axes[0,1].hist(ratios,bins=50,weights=weight,histtype='step',label=variant)
            for ax,key,edges in [(axes[1,0],'installed',PT_EDGES),(axes[1,1],'WTAeta',np.linspace(-1.6,1.6,9))]:
                profile=weighted_profile(np.array([r[key] for r in rows]),ratios,weight,edges)
                valid=[q for q in profile if q['mean'] is not None]
                ax.plot([(q['low']+q['high'])/2 for q in valid],[q['mean'] for q in valid],'.-',label=variant)
                ax.set(xlabel=key,ylabel='mean candidate/installed')
        axes[0,0].set(xlabel='jet pT [GeV]',ylabel='common-event weighted jets');axes[0,0].legend(fontsize=7)
        axes[0,1].set(xlabel='candidate/installed',ylabel='weighted jets');axes[0,1].legend(fontsize=7)
        stem=out/(label+'_isolated_residual_and_l2');save_fig(fig,stem);artifacts.append(str(stem)+'.pdf')
    groups={'selected':(-math.inf,math.inf)} if system=='pp' else GROUPS
    variables=['rawpt','stored_jtpt','installed','candidate_raw_input']
    colors=['0.55','tab:green','tab:blue','tab:orange']
    for group,(low,high) in groups.items():
        chosen=[r for r in rows if (system=='pp' or low<=r['hiBin']<high)]
        w=np.array([r['event_weight'] for r in chosen]);prefix=label+'_'+group
        fig,axes=plt.subplots(2,2,figsize=(11,7.5))
        fig.suptitle(prefix+': paired R4 payload diagnostics; candidate NOT approved; n=1,2 share jets',fontsize=10)
        for v,color in zip(variables,colors):
            values=np.array([r[v] for r in chosen]);hist,edges=np.histogram(values,bins=np.linspace(0,250,51),weights=w)
            axes[0,0].step(edges[:-1],hist,where='post',label=v,color=color)
            h=ROOT.TH1D(prefix+'_'+v,v,50,0,250);h.Sumw2()
            for x,weight in zip(values,w):h.Fill(float(x),float(weight))
            h.Write()
        axes[0,0].set(xlabel='jet pT [GeV]',ylabel='sum of common event weights');axes[0,0].legend(fontsize=7)
        ratios=np.array([r['candidate_over_installed'] for r in chosen]);axes[0,1].hist(ratios,bins=50,weights=w,histtype='step')
        axes[0,1].set(xlabel='candidate(raw input)/installed',ylabel='weighted jets')
        lookup=np.array([r['coordinate_only']/r['installed'] for r in chosen]);axes[1,0].hist(lookup,bins=50,weights=w,histtype='step',color='tab:purple')
        axes[1,0].set(xlabel='standard/WTA lookup: SAME installed payload',ylabel='weighted jets');axes[1,0].set_title('Coordinate-only diagnostic; WTA cones fixed',fontsize=9)
        migration=np.zeros((5,5))
        for r in chosen:migration[r['old_bin'],r['candidate_bin']]+=r['event_weight']
        image=axes[1,1].imshow(migration,origin='lower',cmap='viridis');fig.colorbar(image,ax=axes[1,1],shrink=.75)
        labels=['<40','40–80','80–120','120–160','≥160']
        axes[1,1].set_xticks(range(5));axes[1,1].set_xticklabels(labels,rotation=30,fontsize=7)
        axes[1,1].set_yticks(range(5));axes[1,1].set_yticklabels(labels,fontsize=7)
        axes[1,1].set(xlabel='candidate pT bin',ylabel='installed pT bin')
        stem=out/(prefix+'_pt_ratios_migrations');save_fig(fig,stem);artifacts.append(str(stem)+'.pdf')
        profile_group={'jet_count':len(chosen),'migration_weighted':migration.tolist(),'response':{}}
        if not data:
            fig,axes=plt.subplots(3,2,figsize=(11,10));fig.suptitle(prefix+': common matches BEFORE analysis pT cuts; low-statistics cells remain visible',fontsize=10)
            for rownum,(coordinate,edges) in enumerate([('genpt',PT_EDGES[1:]),('geneta',np.linspace(-1.6,1.6,9)),('centrality',np.array([0,10,30,50,70,90,100.]))]):
                if system=='pp' and coordinate=='centrality':
                    for ax in axes[rownum]:ax.text(.5,.5,'pp activity coordinate not established\nNOT_RUN',ha='center',va='center',transform=ax.transAxes)
                    continue
                x=np.array([r[coordinate] for r in chosen]);gen=np.array([r['genpt'] for r in chosen])
                for variable,color in zip(variables,colors):
                    response=np.array([r[variable] for r in chosen])/gen
                    values=weighted_profile(x,response,w,edges)
                    profile_group['response'][coordinate+'_'+variable]=values
                    for col,statistic in enumerate(('mean','width')):
                        valid=[r for r in values if r[statistic] is not None]
                        axes[rownum,col].plot([(r['low']+r['high'])/2 for r in valid],[r[statistic] for r in valid],'.-',label=variable,color=color)
                        axes[rownum,col].set(xlabel=coordinate,ylabel='response '+statistic)
                    hist=ROOT.TH2D(prefix+'_'+coordinate+'_'+variable,'raw response; no independent-jet error claim',len(edges)-1,array('d',edges),100,0,3)
                    for xi,yi,wi in zip(x,response,w):hist.Fill(float(xi),float(yi),float(wi))
                    hist.Write()
            axes[0,0].legend(fontsize=7);stem=out/(prefix+'_response_mean_width');save_fig(fig,stem);artifacts.append(str(stem)+'.pdf')
        profiles[group]=profile_group
    # Actual before map only: no fake zero-loss or unchanged-after curve.
    fig,ax=plt.subplots(figsize=(7,5))
    if rows:
        graph=ax.hist2d([r['jteta'] for r in rows],[r['jtphi'] for r in rows],bins=[32,36],weights=[r['event_weight'] for r in rows],range=[[-1.6,1.6],[-math.pi,math.pi]])
        fig.colorbar(graph[3],ax=ax,label='matched jet event weights')
    ax.set(xlabel='standard detector jteta',ylabel='standard detector jtphi',title=label+': before veto only; applicable map/after NOT_RUN')
    stem=out/(label+'_eta_phi_before_veto');save_fig(fig,stem);artifacts.append(str(stem)+'.pdf')
    rootfile.Close();return profiles,artifacts


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--oo',type=Path,required=True);p.add_argument('--pp',type=Path,required=True)
    p.add_argument('--oo-data',type=Path);p.add_argument('--pp-data',type=Path)
    p.add_argument('--oo-events',type=int,default=10000);p.add_argument('--pp-events',type=int,default=3000)
    p.add_argument('--upstream',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if not all(1<=n<=10000 for n in (args.oo_events,args.pp_events)):p.error('Use 1..10000 events per bounded sample')
    args.output.mkdir(parents=True,exist_ok=False);ROOT.gROOT.SetBatch(True)
    command=setup_library(args.output);coeff=read_weight_coefficients();results=[]
    for system,path,cap,data in [('OO',args.oo,args.oo_events,False),('pp',args.pp,args.pp_events,False),('OO',args.oo_data,min(args.oo_events,1000),True),('pp',args.pp_data,min(args.pp_events,1000),True)]:
        if path is None:continue
        try:
            rows,metadata=scan(system,path,cap,args.output,args.upstream,coeff,data)
        except ValueError as error:
            if not data:
                raise
            results.append({'label':system+'_data','system':system,'data':True,
                            'status':'BLOCKED_INPUT_SCHEMA','input':str(path.resolve()),
                            'input_sha256':digest(path),'reason':str(error),'plots':[],
                            'comparison':'NOT_RUN; do not substitute jtptUncorrected or another axis without producer equivalence'})
            print(system+'_data BLOCKED_INPUT_SCHEMA: '+str(error),flush=True)
            continue
        profiles,plots=make_plots(rows,metadata,args.output)
        metadata.update(profiles=profiles,plots=plots);results.append(metadata)
        print(metadata['label'],len(rows),'common jets; candidate NOT_APPROVED',flush=True)
    result={'status':'BOUNDED_PAIRED_DIAGNOSTICS_NOT_PHYSICS_CLOSURE','samples':results,'weight_coefficients_from_current_source':coeff,
            'compile_command':command,'head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            'source_sha256':{str(x.relative_to(REPO)):digest(x) for x in [Path(__file__).resolve(),REPO/'src/JetCorrector.cxx',REPO/'src/JetCorrector.h',REPO/'src/EECAnalyzer.cxx',REPO/'src/HighForestReader.cxx']},
            'unrun':{'standalone_L1_plus_V2':'BLOCKED producer/rho/area applicability; no fabricated rho',
                     'data_residual_response':'not an MC correction; no gen truth in data',
                     'event_veto_loss_after_maps':'BLOCKED era/run/collection/map-key/ID/pT prescription approval',
                     'EEC_n1_n2_closure':'NOT_RUN; these jet-level diagnostics share the same input jets and event weights for both powers',
                     'R8':'raw-only independent workstream; not calibrated by this R4 tool'}}
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
