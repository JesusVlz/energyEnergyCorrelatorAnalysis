#!/usr/bin/env python3
"""Draw checked native responses and executed scans; never fabricate missing plots."""
import argparse
import json
from pathlib import Path
import numpy as np
from study_math import REPORT_BINS, interval_index


def label(record):
    c=record['channel']; lo,hi=c['activity']
    if c['system']=='pp': activity='inclusive ppRef'
    elif lo<0 and hi==100: activity='OO selected-MB inclusive (stored -1..100)'
    else: activity=f'OO {lo:g}-{hi:g}%'
    prefix='SYNTHETIC SOFTWARE TEST; ' if record.get('synthetic') else ''
    return prefix+f'{activity}; R=0.4; particle pT>2 GeV; n={c["power"]}'


def draw(folder,formats=('png','pdf')):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    record=json.loads((folder/'manifest.json').read_text())
    if record.get('status') not in ('PREFLIGHT_ONLY','SCAN_EXECUTED_DIAGNOSTIC_ONLY'):
        raise ValueError('No checked products or executed scan: '+str(folder))
    figures=folder/'figures'; figures.mkdir(exist_ok=False)
    data=np.load(folder/'plot_data.npz',allow_pickle=False)
    title=label(record); jets=data['jets']; dr=data['dr']; n=len(dr)-1
    def save(fig,name):
        for ext in formats: fig.savefig(figures/(name+'.'+ext),dpi=160,bbox_inches='tight')
        plt.close(fig)
    try:
        for s in (0,1):
            c=data[f'response_{s}']; positive=c[c>0]
            if not len(positive): raise ValueError('No positive response entries')
            fig,ax=plt.subplots(figsize=(10,8))
            vmin,vmax=float(positive.min()),float(positive.max())
            if vmax<=vmin: vmax=vmin*1.01
            mesh=ax.pcolormesh(np.arange(c.shape[0]+1),np.arange(c.shape[1]+1),
                np.ma.masked_less_equal(c.T,0),norm=LogNorm(vmin=vmin,vmax=vmax),shading='flat')
            fig.colorbar(mesh,ax=ax,label='EEC-weighted matched entries')
            for b in range(1,len(jets)-1):
                ax.axvline(b*n,lw=.5,alpha=.4); ax.axhline(b*n,lw=.5,alpha=.4)
            ticks=(np.arange(len(jets)-1)+.5)*n
            bins=[f'{a:g}-{b:g}'+(' [guard]' if a<40 or a>=160 else '') for a,b in zip(jets[:-1],jets[1:])]
            ax.set_xticks(ticks,bins,rotation=60,ha='right'); ax.set_yticks(ticks,bins)
            ax.set_xlabel('Reco flattened (Delta r, jet pT) blocks [GeV]')
            ax.set_ylabel('Gen flattened (Delta r, jet pT) blocks [GeV]')
            ax.set_title(title+f'\nFold {s+1} response; native within-block bin order')
            fig.tight_layout(); save(fig,f'response_fold{s+1}')
        if record['status']=='PREFLIGHT_ONLY': return figures
        k=record['candidate']['iterations']
        fig,ax=plt.subplots(figsize=(9,5.5))
        for s in (0,1): ax.plot(np.arange(1,21),data['chi2'][s],marker='o',label=f'Split {s+1}')
        ax.plot(np.arange(1,21),data['chi2_average'],marker='s',label='Split average')
        ax.axvline(k,ls='--',label=f'Diagnostic minimum k={k}')
        ax.set(xlabel="D'Agostini iterations",ylabel='Summed diagonal unity chi2',xticks=np.arange(1,21,2))
        ax.set_title(title); ax.legend()
        fig.text(.5,.01,'Diagnostic only: unfolded/truth cross covariance is not included; no calibrated p-value.',ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.05,1,1)); save(fig,'iterations')
        for shape in (True,False):
            fig,axes=plt.subplots(1,3,figsize=(15,5))
            for ax,bounds in zip(axes,REPORT_BINS):
                j=interval_index(jets,bounds); sl=slice(j*n,(j+1)*n)
                for s in (0,1):
                    if shape: r,e=data[f'ratio_{s}_{k}_{j}'],data[f'error_{s}_{k}_{j}']
                    else:
                        u,t=data[f'u_{s}_{k}'][sl],data[f'truth_{s}'][sl]
                        vu,vt=np.diag(data[f'cov_{s}_{k}'])[sl],data[f'truth_var_{s}'][sl]
                        r=np.divide(u,t,out=np.full_like(u,np.nan),where=t>0)
                        v=np.divide(vu,t*t,out=np.full_like(u,np.nan),where=t>0)
                        v+=np.divide(u*u*vt,t**4,out=np.full_like(u,np.nan),where=t>0)
                        if np.any(v[np.isfinite(v)] < 0): raise ValueError('Negative ratio variance')
                        e=np.sqrt(v)
                    x=np.sqrt(dr[:-1]*dr[1:]); valid=np.isfinite(r)&np.isfinite(e)&(x>0)
                    ax.errorbar(x[valid],r[valid],yerr=e[valid],fmt='o',ms=3,label=f'Split {s+1}, {k} iterations')
                ax.axhline(1,ls='--',lw=1); ax.set_xscale('log'); ax.set_xlabel('Delta r')
                ax.set_title(f'{bounds[0]:g}<jet pT<{bounds[1]:g} GeV')
                ax.set_ylabel('Unfolded / truth'+(' (normalized shape)' if shape else ' (weighted yield)'))
                ax.legend(fontsize=8)
            fig.suptitle(title+'\nTruth: '+record['truth_definition'])
            fig.text(.5,.01,'PILOT: truth errors use a Sumw2 proxy; joint closure covariance pending. No y-range clipping.',ha='center',fontsize=9)
            fig.tight_layout(rect=(0,.06,1,.88)); save(fig,'closure_shape' if shape else 'closure_yield')
        return figures
    finally: data.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path,help='One channel folder or top-level run folder')
    args=parser.parse_args()
    folders=[args.directory] if (args.directory/'plot_data.npz').is_file() else sorted(p.parent for p in args.directory.glob('*/plot_data.npz'))
    if not folders: parser.error('No checked archive; no placeholder plots produced')
    for p in folders: print(draw(p))


if __name__=='__main__': main()
