#!/usr/bin/env python3
"""Explain unresolved positive references in the same bounded input-QA entries."""
import argparse,json,math,sys
from collections import Counter
from pathlib import Path
from audit_lowpt_inputs import REPO,digest,select_branches

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);a=p.parse_args();qa=json.loads((a.directory/'input_qa.json').read_text())
 import ROOT
 ROOT.gROOT.SetBatch(True);ROOT.TH1.AddDirectory(False);ROOT.gInterpreter.AddIncludePath(str(REPO/'src'));ROOT.gInterpreter.Declare('#include "'+str(REPO/'src/JetCorrector.cxx')+'"\n#include "'+str(REPO/'src/Run3JetSelection.cxx')+'"')
 baseline=REPO/'validation/run3Pairs/results/2026-11-02-an751-v1/expanded';tasks=json.loads((baseline/'sample_manifest.json').read_text())['tasks'];reports=[]
 for item in qa['files']:
  source=Path(item['path']).resolve();task=next(t for t in tasks if Path(t['input']).resolve()==source);meta_path=Path(task['output']+'.metadata.json');meta=json.loads(meta_path.read_text());ledger=Path(meta['native_outputs'])/'reco.jsonl.events.jsonl';limit=item['entry_range'][1]
  event_status={}
  with ledger.open() as stream:
   for line in stream:
    e=json.loads(line)
    if e['source_entry']>=limit:break
    event_status[e['source_entry']]=e
  if len(event_status)!=limit:raise ValueError('Incomplete source-event ledger')
  components={};selected_entries=set();native_jets=Counter();product_mismatches=0
  with Path(task['output']).open() as stream:
   for line in stream:
    row=json.loads(line)
    if row['source_entry']>=limit:break
    selected_entries.add(row['source_entry'])
    for key in ['generator_weight','stored_generator_weight','vz_weight','centrality_weight','event_weight']:
     value=float(row[key]);stat=components.setdefault(key,{'count':0,'negative':0,'zero':0,'nonfinite':0,'unity':0,'min':None,'max':None});stat['count']+=1;stat['negative']+=value<0;stat['zero']+=value==0;stat['nonfinite']+=not math.isfinite(value);stat['unity']+=value==1
     if math.isfinite(value):stat['min']=value if stat['min'] is None else min(stat['min'],value);stat['max']=value if stat['max'] is None else max(stat['max'],value)
    product_mismatches+=not math.isclose(row['event_weight'],row['generator_weight']*row['vz_weight']*row['centrality_weight'],rel_tol=1e-12)
    for jet in row['jets']:
     if 40<=jet['pt']<80:native_jets['selected_reco40_80']+=1
  for entry,event in event_status.items():event['event_filters_pass_before_veto']=event['event_filters_pass'];event['event_filters_pass']=entry in selected_entries
  if len(selected_entries)!=item['counts']['selected_events']:raise ValueError('Native selected-event population disagrees with current QA selection')
  f=ROOT.TFile.Open(str(source));j=f.Get(item['jet_tree']);cols=['nref','ngen','rawpt','jtpt','jteta','jtphi','WTAeta','WTAphi','refpt','refeta','refphi','genpt','geneta','genphi','WTAgeneta','WTAgenphi','jtPfCHF','jtPfNHF','jtPfCEF','jtPfNEF','jtPfMUF','jtPfCHM','jtPfNHM','jtPfCEM','jtPfNEM','jtPfMUM'];select_branches(j,cols)
  corrector=ROOT.JetCorrector(item['jec']['path']);prescription=ROOT.EECJetSelection.JetIDPrescription.HINPbPb if item['system']=='OO' else ROOT.EECJetSelection.JetIDPrescription.HINpp
  unresolved=[];crossings=[];tot=Counter()
  for entry in range(limit):
   j.GetEntry(entry);event=event_status[entry]
   for r in range(int(j.nref)):
    refpt=float(j.refpt[r]);matches=[k for k in range(int(j.ngen)) if abs(float(j.refeta[r])-float(j.geneta[k]))<.015 and abs(float(j.refphi[r])-float(j.genphi[k]))<.015 and abs(refpt-float(j.genpt[k]))<.03*float(j.genpt[k])]
    missing=refpt>=0 and not matches
    if not missing and not(event['event_filters_pass'] and abs(float(j.WTAeta[r]))<1.6 and matches):continue
    corrector.SetJetPT(float(j.rawpt[r]));corrector.SetJetEta(float(j.WTAeta[r]));corrector.SetJetPhi(float(j.WTAphi[r]));pt=float(corrector.GetCorrectedPT())
    selected=False
    if event['event_filters_pass'] and abs(float(j.WTAeta[r]))<1.6:
     x=ROOT.EECJetSelection.JetInputs();x.eta=float(j.jteta[r]);x.phi=float(j.jtphi[r]);x.rawPt=float(j.rawpt[r]);x.storedPt=float(j.jtpt[r]);x.neutralHadronFraction=float(j.jtPfNHF[r]);x.neutralEmFraction=float(j.jtPfNEF[r]);x.chargedHadronFraction=float(j.jtPfCHF[r]);x.chargedEmFraction=float(j.jtPfCEF[r]);x.muonFraction=float(j.jtPfMUF[r]);x.chargedMultiplicity=int(j.jtPfCHM[r])+int(j.jtPfCEM[r])+int(j.jtPfMUM[r]);x.neutralMultiplicity=int(j.jtPfNHM[r])+int(j.jtPfNEM[r]);x.numberOfConstituents=x.chargedMultiplicity+x.neutralMultiplicity;x.multiplicitiesAreLowerBounds=True;selected=bool(ROOT.EECJetSelection.PassJetID(x,prescription))
    row={'source_entry':entry,'event_id':event['id'],'reco_index':r,'stored_hiBin':event['activity'],'refpt':refpt,'refeta':float(j.refeta[r]),'refphi':float(j.refphi[r]),'recomputed_reco_pt':pt,'reco_wta_eta':float(j.WTAeta[r]),'selected_event_jet_without_pt_cut':selected,'selected_reco_ge40':selected and pt>=40,'retained_ngen':int(j.ngen),'retained_gen_min_pt':min([float(j.genpt[k]) for k in range(int(j.ngen))],default=None)}
    if missing:
     # This diagnoses identity/reference retention only. Do not replace matching.
     row['reason']='REFERENCE_BELOW_OBSERVED_GEN5_BOUNDARY' if refpt<5 else 'REFERENCE_NOT_RESOLVED_WITHIN_READER_TOLERANCES'
     unresolved.append(row);tot['positive_ref_unresolved']+=1;tot['unresolved_selected_reco_ge40']+=selected and pt>=40;tot['unresolved_ref_below5']+=refpt<5
    elif selected and pt>=40:
     k=matches[0];gp=float(j.genpt[k]);eta=float(j.WTAgeneta[k]);tot['selected_reco_ge40_matched']+=1
     if abs(eta)>=1.6 or gp<10 or gp>=5020:
      row.update({'gen_index':k,'gen_pt':gp,'gen_wta_eta':eta,'outside_truth_eta':abs(eta)>=1.6,'outside_truth_pt_10_5020':not 10<=gp<5020});crossings.append(row);tot['selected_reco_ge40_truth_eta_crossing']+=abs(eta)>=1.6;tot['selected_reco_ge40_truth_pt_crossing']+=not 10<=gp<5020
  if len(unresolved)!=item['counts']['positive_reference_unresolved']:raise ValueError('Targeted audit count differs')
  reports.append({'system':item['system'],'path':str(source),'entry_range':item['entry_range'],'native_event_ledger':str(ledger),'native_metadata_sha256':digest(meta_path),'event_selection_crosscheck':'PASS exact selected-event counts from merged native export against current raw QA; scanned ledger flag alone precedes veto','native_merged_export':task['output'],'native_event_weight_components':components,'native_event_weight_product_mismatches':product_mismatches,'native_jet_counts':native_jets,'raw_activity_range':[min(e['activity'] for e in event_status.values()),max(e['activity'] for e in event_status.values())],'raw_invalid_activity_count':sum(not 0<=e['activity']<200 for e in event_status.values()) if item['system']=='OO' else None,'counts':tot,'unresolved':unresolved,'reco_ge40_truth_acceptance_crossings':crossings});f.Close()
 result={'status':'FAILURE_SPECIFIC_INPUT_QA_NOT_RESPONSE','weight_audit':'Existing analyzer-native exports for exactly the same source entries; stored generator vs applied generator/vz/activity/event factors separately counted. No component refit. Applied Run3 generator unity is explicit current policy, not a cross-section normalization validation.','purpose':'Explain positive references unresolved in the same20k bounded entries; audit acceptance crossings for selected reco>=40. No matching behavior changed.','source_sha256':digest(Path(__file__)),'input_qa_sha256':digest(a.directory/'input_qa.json'),'files':reports}
 (a.directory/'matching_reference_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({s:dict(sum((r['counts'] for r in reports if r['system']==s),Counter())) for s in ['OO','pp']}))
if __name__=='__main__':main()
