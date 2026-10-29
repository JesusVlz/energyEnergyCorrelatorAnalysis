#!/usr/bin/env python3
"""Compare actual capped ordinary OO fills against disjoint centrality binning.
Run standard manager projection/reload on OO and pp; no physics closure claim.
"""
import argparse
import json
import math
from pathlib import Path
import ROOT
from check_projection import configure


def sparse_values(h):
    return [(h.GetBinContent(i),h.GetBinError(i)) for i in range(h.GetNbins())]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    ROOT.gROOT.SetBatch(True)
    repo=Path(__file__).resolve().parents[2]
    assert ROOT.gSystem.Load(str(args.directory/'projection-fixtures/standardManager.so'))>=0
    ROOT.gInterpreter.AddIncludePath(str(repo))
    ROOT.gInterpreter.Declare('#include "plotting/EECHistogramManager.h"')
    reports={}
    sources={}
    for system in ('OO','pp','OO_split'):
        path=args.directory/(system+'_R04_standard.root')
        if not path.exists():continue
        f=ROOT.TFile(str(path));sources[system]=f
        h=f.Get('energyEnergyCorrelator');q=f.Get('eventQA/eventCoverage')
        assert h and q and not f.Get('inclusive')
        reports[system]={'status':'BOUNDED_STANDARD_EEC_SOFTWARE_DIAGNOSTIC',
                         'selected_events':sum(q.GetBinContent(x,4) for x in range(1,6)),
                         'read_events':sum(q.GetBinContent(x,1) for x in range(1,6)),
                         'pair_entries':h.GetEntries(), 'pairs_sumw':sum(v for v,e in sparse_values(h)),
                         'pairs_sumw2':sum(e*e for v,e in sparse_values(h))}
        assert reports[system]['read_events']==100
        if system=='OO_split':continue
        reports[system]['projected']={}
        for n in (1,2):
            m=configure(f,n);m.LoadHistograms()
            expected=[]
            for j in range(3):
                hist=m.GetHistogramEnergyEnergyCorrelator(0,0,j,3,0,4)
                expected.append(hist.Integral('width'))
            output=args.directory/f'{system}_R04_N{n}_projected.root'
            m.Write(str(output),'RECREATE')
            loaded=ROOT.TFile(str(output));restored=configure(loaded,n,True)
            restored.LoadProcessedHistograms()
            for j,total in enumerate(expected):
                assert math.isclose(restored.GetHistogramEnergyEnergyCorrelator(0,0,j,3,0,4).Integral('width'),total,rel_tol=1e-7,abs_tol=1e-10)
            reports[system]['projected'][str(n)]=expected
            del restored
            loaded.Close()
            del m
    if 'OO_split' in reports:
        for key in ('selected_events','read_events','pair_entries','pairs_sumw','pairs_sumw2'):
            assert math.isclose(reports['OO'][key],reports['OO_split'][key],rel_tol=1e-7,abs_tol=1e-9),(key,reports['OO'][key],reports['OO_split'][key])
        for name in ('eventQA/eventActivity',):
            assert sparse_values(sources['OO'].Get(name))==sparse_values(sources['OO_split'].Get(name)),name
        reports['binning_invariance']='PASS: identical counts/activity ledger; summed pair yields and Sumw2 agree within 1e-7 relative Float_t accumulation tolerance over same accepted range'
    args.output.write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))
    for f in sources.values():f.Close()


if __name__=='__main__':main()
