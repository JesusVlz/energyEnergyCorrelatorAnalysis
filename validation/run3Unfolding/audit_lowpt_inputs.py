#!/usr/bin/env python3
"""Read-only bounded raw-input QA and embedded-metadata inventory; no response filling."""
import argparse,csv,hashlib,json,math,sys,time
from array import array
from collections import Counter,defaultdict
from pathlib import Path
REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'validation/run3Pairs'))
from export_forest_legacy import select_branches,digest
PT=[0,10,20,40,80,120,160,200,300,500,1000,5020]

def card(f):
 d=f.Get('JCard');out={}
 if not d:return out
 for name in sorted({k.GetName() for k in d.GetListOfKeys()}):
  x=d.Get(name)
  if x.InheritsFrom('TVectorT<float>'):out[name]=[float(x.GetMatrixArray()[i]) for i in range(x.GetNrows())]
  elif x.InheritsFrom('TObjString'):out[name]=str(x.GetString().Data())
 return out

def axis(a):return {'title':str(a.GetTitle()),'edges':[float(a.GetBinLowEdge(i)) for i in range(1,a.GetNbins()+2)]}
def objects(d,prefix=''):
 for name in sorted({k.GetName() for k in d.GetListOfKeys()}):
  if not prefix and not any(s in name.lower() for s in ('unfold','response','eventqa','inclusivejet','energyenergycorrelator')):continue
  x=d.Get(name);p=prefix+name
  if x.InheritsFrom('TDirectory'):
   yield from objects(x,p+'/')
  elif x.InheritsFrom('THnSparse'):
   yield {'path':p,'class':x.ClassName(),'entries':float(x.GetEntries()),'axes':[axis(x.GetAxis(i)) for i in range(x.GetNdimensions())]}
  elif x.InheritsFrom('TH1'):
   axes=[x.GetXaxis(),x.GetYaxis(),x.GetZaxis()][:x.GetDimension()]
   yield {'path':p,'class':x.ClassName(),'entries':float(x.GetEntries()),'axes':[axis(a) for a in axes]}

def inventory(ROOT,manifest,out):
 rows=[];details={'raw_manifest':str(manifest),'raw_manifest_sha256':digest(manifest),'products':[],'axis_schemas':{}}
 m=json.loads(manifest.read_text())
 for t in m['tasks']:
  p=Path(t['input']);s=p.stat();frozen=t['input_identity']
  rows.append({'kind':'raw_MC_forest','path':str(p),'system':t['system'],'events':t['entry_limit'],'metadata_source':'frozen header manifest; schema QA only representative files','status':'PRESENT_FROZEN_IDENTITY' if s.st_size==frozen['size_bytes'] and s.st_mtime_ns==frozen['mtime_ns'] else 'CHANGED','size_bytes':s.st_size,'schema_or_card':'input_qa.json for bounded representative files; upstream retention completeness unresolved'})
 paths=set(REPO.glob('*.root'))|set((REPO/'projected').glob('*.root'))
 for d in ['figures/OO_inclusive_Run3_20261029','figures/OO_inclusive_lowPt_Run3_20261029','validation/run3Jets/results/2026-10-29-mc-veto']:
  paths.update((REPO/d).glob('*.root'))
 paths.update((REPO/'validation/run3Pairs/results/2026-11-02-an751-v1/expanded').glob('*.native/*.root'))
 for p in sorted(paths):
  f=ROOT.TFile.Open(str(p));c=card(f);h=list(objects(f));nonempty=[x['path'] for x in h if 'unfold' in x['path'].lower() and x['entries']>0]
  val=lambda k:c.get(k,[None])[0] if isinstance(c.get(k,[None]),list) else c[k]
  system={9:'OO_data',10:'OO_MC',11:'ppRef_data',12:'ppRef_MC'}.get(val('DataType'),'unknown')
  row={'kind':'processed_ROOT','path':str(p),'system':system,'size_bytes':p.stat().st_size,'metadata_source':'embedded JCard and actual histogram objects','events':None,'status':'AUDITED_NOT_CERTIFIED_RESPONSE','datatype':val('DataType'),'mc_mode':val('McCorrelationType'),'radius':val('JetRadius'),'jet_axis':val('JetAxis'),'min_reco_unfold':val('MinJetPtUnfoldingReco'),'min_truth_unfold':val('MinJetPtUnfoldingTruth'),'min_pair_track_pt':val('MinPairTrackPt'),'weight_exponents':json.dumps(c.get('WeightExponent')),'filled_histograms':val('FilledHistograms'),'run3_jet_id':val('Run3JetID'),'jet_veto':val('ApplyJetVeto'),'trigger':val('TriggerSelection'),'hf_shift':val('HFEnergyShift'),'nonempty_unfold_objects':len(nonempty),'schema_or_card':json.dumps({k:c.get(k) for k in ['CentralityBinEdges','JetPtBinEdgesEEC','TrackPtBinEdgesEEC','InputFile','GitHash']})}
  n=f.Get('nEvents')
  if n and n.InheritsFrom('TH1'):row['events']=float(n.GetBinContent(1))
  for item in h:
   payload=json.dumps(item.pop('axes'),sort_keys=True);key=hashlib.sha256(payload.encode()).hexdigest()
   details['axis_schemas'].setdefault(key,json.loads(payload));item['axis_schema_sha256']=key
  details['products'].append({'path':str(p),'card':c,'histograms':h,'nonempty_unfold_objects':nonempty})
  rows.append(row);f.Close()
 fields=sorted(set().union(*(r.keys() for r in rows)))
 with (out/'input_inventory.csv').open('w') as stream:
  w=csv.DictWriter(stream,fields);w.writeheader();w.writerows(rows)
 (out/'input_inventory_details.json').write_text(json.dumps(details,indent=2,allow_nan=False)+'\n')
 return {'raw_files':len(m['tasks']),'raw_entries':dict(Counter({s:sum(t['entry_limit'] for t in m['tasks'] if t['system']==s) for s in ['OO','pp']})),'products':len(paths)}

def string_leaf(t,name):
 leaf=t.GetLeaf(name)
 try:return str(leaf.GetValueString())
 except Exception:return 'UNRESOLVED_STRING_LEAF'

def channels(system,h):
 if system=='pp':return ['ppRef']
 if not 0<=h<200:return ['selected_MB','invalid_activity']
 return ['selected_MB', 'C0_30' if h<60 else 'C30_50' if h<100 else 'C50_90' if h<180 else 'C90_100']

def qa(ROOT,tasks,out,limit,files_per_system):
 ROOT.gInterpreter.AddIncludePath(str(REPO/'src'))
 ROOT.gInterpreter.Declare('#include "'+str(REPO/'src/JetCorrector.cxx')+'"\n#include "'+str(REPO/'src/Run3JetSelection.cxx')+'"')
 summaries=[];hist={};root=ROOT.TFile(str(out/'input_qa.root'),'RECREATE');seen=set()
 def geth(name,title,bins=PT):
  if name not in hist:
   hist[name]=ROOT.TH1D(name,title,len(bins)-1,array('d',bins));hist[name].SetDirectory(0)
  return hist[name]
 for system in ['OO','pp']:
  jec_path=REPO/'jetEnergyCorrections'/('L2Relative_AK4PF_OO_V1.txt' if system=='OO' else 'Prompt24HIpp_V1_MC_L2Relative_AK4PF.txt')
  veto_path=REPO/'jetVetoMaps'/('OO2025/Summer24Prompt25_RunCDEFG.root' if system=='OO' else 'ppRef2024/Summer24Prompt24_RunBCDEFGHI.root')
  corrector=ROOT.JetCorrector(str(jec_path));veto=ROOT.EECJetSelection.JetVetoMap(str(veto_path),'jetvetomap_all');prescription=ROOT.EECJetSelection.JetIDPrescription.HINPbPb if system=='OO' else ROOT.EECJetSelection.JetIDPrescription.HINpp
  for task in [t for t in tasks if t['system']==system][:files_per_system]:
   p=Path(task['input']);f=ROOT.TFile.Open(str(p));e=f.Get('hiEvtAnalyzer/HiTree');jetname=('akCs4' if system=='OO' else 'ak4')+'PFJetAnalyzer/t';j=f.Get(jetname);skim=f.Get('skimanalysis/HltTree');g=f.Get('HiGenParticleAna/hi');trk=f.Get('PbPbTracks/trackTree' if system=='OO' else 'ppTracks/trackTree');info=f.Get('HiForestInfo/HiForest');info.GetEntry(0)
   count=min(limit,int(e.GetEntries()));c=Counter();ch=defaultdict(Counter);ranges={};examples=[]
   summary={'system':system,'path':str(p),'entry_range':[0,count],'full_file_entries':int(e.GetEntries()),'sample_namespace':task['sample_namespace'],'producer_info':{b.GetName():string_leaf(info,b.GetName()) for b in info.GetListOfBranches()},'native_jet_collections':[k.GetName() for k in f.GetListOfKeys() if 'JetAnalyzer' in k.GetName()],'jet_tree':jetname,'jet_branches':[b.GetName() for b in j.GetListOfBranches()],'gen_particle_branches':[b.GetName() for b in g.GetListOfBranches()],'track_branches':[b.GetName() for b in trk.GetListOfBranches()],'tree_entries':{name:int(f.Get(name).GetEntries()) for name in ['hiEvtAnalyzer/HiTree',jetname,'skimanalysis/HltTree','HiGenParticleAna/hi','PbPbTracks/trackTree' if system=='OO' else 'ppTracks/trackTree']},'jec':{'path':str(jec_path),'sha256':digest(jec_path)},'veto':{'path':str(veto_path),'sha256':digest(veto_path),'histogram':'jetvetomap_all','pt_convention':'stored'},'counts':c,'channels':ch,'observed_ranges':ranges,'examples':examples}
   jetcols=['nref','ngen','rawpt','jtpt','jteta','jtphi','WTAeta','WTAphi','refpt','refeta','refphi','refdrjt','genpt','geneta','genphi','WTAgeneta','WTAgenphi','genmatchindex','jtPfCHF','jtPfNHF','jtPfCEF','jtPfNEF','jtPfMUF','jtPfCHM','jtPfNHM','jtPfCEM','jtPfNEM','jtPfMUM']
   select_branches(j,jetcols);select_branches(e,['run','lumi','evt','vz','weight','hiBin','pthat']);filters=['pprimaryVertexFilter']+(['OOphfCoincFilterPF2Th4','pileupVertexFilter'] if system=='OO' else []);select_branches(skim,filters)
   def rangefill(name,v):
    if math.isfinite(v):
     old=ranges.setdefault(name,{'min':v,'max':v});old['min']=min(old['min'],v);old['max']=max(old['max'],v)
    else:c[name+'_nonfinite']+=1
   for n in range(count):
    e.GetEntry(n);j.GetEntry(n);skim.GetEntry(n);c['raw_events']+=1
    identity=(task['sample_namespace'],int(e.run),int(e.lumi),int(e.evt));c['duplicate_event_ids']+=identity in seen;seen.add(identity)
    w=float(e.weight);c['stored_weight_negative']+=w<0;c['stored_weight_zero']+=w==0;c['stored_weight_nonfinite']+=not math.isfinite(w);rangefill('stored_weight',w);rangefill('pthat',float(e.pthat))
    h=int(e.hiBin);rangefill('stored_hiBin',h);c['invalid_hiBin']+=system=='OO' and not 0<=h<200;chs=channels(system,h)
    inputs=[];passing=[];vetoed=False
    for r in range(int(j.nref)):
     x=ROOT.EECJetSelection.JetInputs();x.eta=float(j.jteta[r]);x.phi=float(j.jtphi[r]);x.rawPt=float(j.rawpt[r]);x.storedPt=float(j.jtpt[r]);x.neutralHadronFraction=float(j.jtPfNHF[r]);x.neutralEmFraction=float(j.jtPfNEF[r]);x.chargedHadronFraction=float(j.jtPfCHF[r]);x.chargedEmFraction=float(j.jtPfCEF[r]);x.muonFraction=float(j.jtPfMUF[r]);x.chargedMultiplicity=int(j.jtPfCHM[r])+int(j.jtPfCEM[r])+int(j.jtPfMUM[r]);x.neutralMultiplicity=int(j.jtPfNHM[r])+int(j.jtPfNEM[r]);x.numberOfConstituents=x.chargedMultiplicity+x.neutralMultiplicity;x.multiplicitiesAreLowerBounds=True
     inputs.append(x);passing.append(bool(ROOT.EECJetSelection.PassJetID(x,prescription)))
     vetoed=vetoed or bool(ROOT.EECJetSelection.RejectEventForJet(x,x.storedPt,veto,prescription))
     rangefill('stored_reco_pt',x.storedPt);rangefill('raw_reco_pt',x.rawPt)
    ev=all(bool(getattr(skim,b)) for b in filters) and abs(float(e.vz))<=15 and not vetoed
    c['event_vetoed']+=vetoed;c['selected_events']+=ev
    for label in chs:ch[label]['raw_events']+=1;ch[label]['selected_events']+=ev
    matches=[];selectedmatches=[]
    for r,x in enumerate(inputs):
     c['stored_reco_jets']+=1;rp=float(j.refpt[r]);c['positive_reference']+=rp>=0
     candidates=[];wrapped=[]
     for k in range(int(j.ngen)):
      if abs(float(j.refeta[r])-float(j.geneta[k]))<.015 and abs(rp-float(j.genpt[k]))<.03*float(j.genpt[k]):
       delta=float(j.refphi[r])-float(j.genphi[k])
       if abs(delta)<.015:candidates.append(k)
       if abs(math.remainder(delta,2*math.pi))<.015:wrapped.append(k)
     k=candidates[0] if candidates else -1;matches.append(k)
     c['unmatched_reco_reference_lookup']+=k<0;c['positive_reference_unresolved']+=rp>=0 and k<0;c['ambiguous_reference_candidates']+=len(candidates)>1;c['phi_wrap_only_matches']+=not candidates and bool(wrapped)
     if len(candidates)>1 and len(examples)<10:examples.append({'entry':n,'reco_index':r,'candidates':candidates})
     corrector.SetJetPT(x.rawPt);corrector.SetJetEta(float(j.WTAeta[r]));corrector.SetJetPhi(float(j.WTAphi[r]));pt=float(corrector.GetCorrectedPT());rangefill('recomputed_reco_pt',pt)
     selected=ev and passing[r] and abs(float(j.WTAeta[r]))<1.6
     if not selected:continue
     c['selected_reco_no_pt_cut']+=1;selectedmatches.append(k)
     if k>=0:
      c['selected_reco_matched_gen_outside_eta']+=abs(float(j.WTAgeneta[k]))>=1.6
      c['selected_reco_matched_gen_below10']+=float(j.genpt[k])<10
      c['selected_reco_matched_gen_above5020']+=float(j.genpt[k])>=5020
      dr=math.hypot(float(j.jteta[r])-float(j.geneta[k]),math.remainder(float(j.jtphi[r])-float(j.genphi[k]),2*math.pi));geth(system+'_matched_dR','Reference-linked selected jets;physical #DeltaR(reco,gen);Jets',[i*.01 for i in range(101)]).Fill(dr);c['refdrjt_mismatch_gt1e-4']+=abs(dr-float(j.refdrjt[r]))>1e-4
     for label in chs:
      ch[label]['reco_no_pt_cut']+=1;ch[label]['reco_matched' if k>=0 else 'reco_unmatched']+=1
      for lo,hi in zip(PT,PT[1:]):
       if lo<=pt<hi:ch[label][f'reco_{lo}_{hi}']+=1
      geth(system+'_'+label+'_reco','Selected events, PF Jet ID, |WTA #eta|<1.6;recomputed jet p_{T} (GeV);Jets').Fill(pt)
    c['duplicate_reference_gen_assignments']+=sum(v-1 for k,v in Counter(matches).items() if k>=0 and v>1)
    c['duplicate_selected_gen_assignments']+=sum(v-1 for k,v in Counter(selectedmatches).items() if k>=0 and v>1)
    for k in range(int(j.ngen)):
     c['stored_gen_jets']+=1;gp=float(j.genpt[k]);rangefill('stored_gen_pt',gp);idx=int(j.genmatchindex[k]);c['genmatchindex_invalid']+=idx<-1 or idx>=int(j.nref)
     c['genmatchindex_disagrees_reference']+=idx>=0 and idx<int(j.nref) and matches[idx]!=k
     c['gen_unmatched_reference_lookup']+=k not in matches
     if not(ev and abs(float(j.WTAgeneta[k]))<1.6):continue
     for label in chs:
      ch[label]['gen_no_pt_cut']+=1
      for lo,hi in zip(PT,PT[1:]):
       if lo<=gp<hi:ch[label][f'gen_{lo}_{hi}']+=1
      geth(system+'_'+label+'_truth','Selected events, |WTA #eta|<1.6;generator jet p_{T} (GeV);Jets').Fill(gp)
   # Read only vector-level retention summaries, not particle-pair responses.
   for name,tree,ptcol,etacol in [('reco_tracks',trk,'trkPt','trkEta'),('gen_particles',g,'pt','eta')]:
    df=ROOT.RDataFrame(tree).Range(count).Define('qa_min',f'{ptcol}.empty() ? -1.f : ROOT::VecOps::Min({ptcol})').Define('qa_max',f'{ptcol}.empty() ? -1.f : ROOT::VecOps::Max({ptcol})').Define('qa_below2',f'ROOT::VecOps::Sum({ptcol}<2.f)').Define('qa_size',f'(unsigned long long){ptcol}.size()').Define('qa_eta_max',f'{etacol}.empty() ? 0.f : ROOT::VecOps::Max(abs({etacol}))')
    vals=df.AsNumpy(['qa_min','qa_max','qa_below2','qa_size','qa_eta_max'])
    summary[name]={'observed_min_pt':min(float(v) for v in vals['qa_min'] if v>=0),'observed_max_pt':max(float(v) for v in vals['qa_max']),'stored_below2':sum(int(v) for v in vals['qa_below2']),'stored_count':sum(int(v) for v in vals['qa_size']),'observed_max_abs_eta':max(float(v) for v in vals['qa_eta_max']),'retention_completeness':'UNPROVEN: observation does not establish producer thresholds or signal-subevent retention'}
   summaries.append(summary);f.Close();print(json.dumps({'system':system,'file':str(p),'counts':dict(c)}),flush=True)
 report={'status':'BOUNDED_INPUT_QA_ONLY_NOT_RESPONSE_OR_CLOSURE','selection':'No centrality cut for inclusive OO or pp; activity partitions are bookkeeping only and invalid activity remains in inclusive counts. TriggerSelection0; native MC event filters and |vz|<=15; exact existing HIN PF Jet ID/event veto helpers on detector coordinates; stored-pT veto; nominal existing L2Relative JEC evaluated on WTA axes as analyzer currently does; no JER (JetUncertainty0). Histogram populations are unweighted and have no pT cut. No pair correction or response matrix is filled.','matching':'Exact current UnfoldingForestReader ref identity: first gen with |deta|<.015, unwrapped |dphi|<.015, |dpT|<.03 genpt; wrapped variant is QA only, not a changed mapping. Producer reco-gen association criterion remains separate.','raw_weight_policy':'Stored generator weights inspected only; current native Run3 generator factor is unity; these QA histograms are unweighted counts. No fitted vz/activity weights applied.','guard_edges':PT,'files':summaries,'histogram_flows':{k:{'underflow':v.GetBinContent(0),'overflow':v.GetBinContent(v.GetNbinsX()+1),'entries':v.GetEntries()} for k,v in hist.items()},'source_sha256':digest(Path(__file__))}
 root.cd()
 for h in hist.values():h.Write()
 ROOT.TObjString(json.dumps(report,allow_nan=False)).Write('InputQAProvenance');root.Close()
 (out/'input_qa.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 (out/'figures').mkdir(exist_ok=True)
 for system in ['OO','pp']:
  names=[c for c in ['selected_MB','C0_30','C30_50','C50_90','C90_100'] if system+'_'+c+'_reco' in hist] if system=='OO' else ['ppRef']
  fig,axs=plt.subplots(1,len(names),figsize=(max(6,3.8*len(names)),4),squeeze=False)
  for ax,label in zip(axs[0],names):
   for pop,color in [('reco','tab:blue'),('truth','tab:orange')]:
    h=hist.get(system+'_'+label+'_'+pop)
    if h:ax.stairs([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)],PT,label=pop,color=color)
   ax.set(xscale='log',yscale='log',xlim=(5,1000),xlabel='Jet pT (GeV)',ylabel='Unweighted jets',title=label);ax.legend();ax.grid(alpha=.2)
  fig.suptitle(system+' bounded input QA: selected events, |WTA eta|<1.6; not a response');fig.tight_layout();fig.savefig(out/'figures'/('input_qa_'+system+'_jet_support.png'),dpi=160);plt.close(fig)
 return report

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--manifest',type=Path,default=REPO/'validation/run3Pairs/results/2026-11-03-an751-all-mc-v1/expanded/sample_manifest.json');p.add_argument('--max-events',type=int,default=5000);p.add_argument('--files-per-system',type=int,default=2);p.add_argument('--inventory-only',action='store_true');p.add_argument('--qa-only',action='store_true');a=p.parse_args()
 import ROOT
 ROOT.gROOT.SetBatch(True);ROOT.TH1.AddDirectory(False);a.output.mkdir(parents=True,exist_ok=True)
 started=time.monotonic();inv=inventory(ROOT,a.manifest,a.output) if not a.qa_only else None
 if not a.inventory_only:qa(ROOT,json.loads(a.manifest.read_text())['tasks'],a.output,a.max_events,a.files_per_system)
 print(json.dumps({'inventory':inv,'elapsed_seconds':time.monotonic()-started,'output':str(a.output.resolve())}))
if __name__=='__main__':main()
