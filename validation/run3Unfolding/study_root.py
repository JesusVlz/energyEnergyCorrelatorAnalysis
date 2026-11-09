"""Thin PyROOT adapter to the existing EEC classes and RooUnfold, imported lazily."""
from array import array
from pathlib import Path
import hashlib
import numpy as np
from study_math import covariance, flattened_layout, interval_index, physical_index

REPO = Path(__file__).resolve().parents[2]
CONTRACT_KEYS = ('DataType','JetRadius','JetType','JetAxis','JetEtaCut','TrackEtaCut',
    'MinTrackPtCut','MaxTrackPtCut','MinPairTrackPt','SubeventCut','Run3JetID','ApplyJetVeto',
    'JetVetoMapFile','JetVetoMapHistogram','JetVetoPtConvention','JetUncertainty',
    'JetPtWeight','SmearDeltaR','SmearEnergyWeight','HFEnergyShift','CentralityBinEdges',
    'JetPtBinEdgesEEC','TrackPtBinEdgesEEC','WeightExponent','BackgroundMethods',
    'ZVertexCut','LowPtHatCut','HighPtHatCut','MultiplicityMode','TriggerSelection',
    'MinMaxTrackPtFraction','MaxMaxTrackPtFraction','JetVetoInputMinPt')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def load_root(drawing=None, library='libRooUnfold'):
    try:
        import ROOT
    except ImportError as error:
        raise RuntimeError('PyROOT unavailable; use the established ROOT/CMSSW runtime') from error
    ROOT.gROOT.SetBatch(True)
    if ROOT.gSystem.Load(library)<0:
        raise RuntimeError('Cannot load RooUnfold: '+library)
    if ROOT.gSystem.Load(str(drawing or REPO/'plotting/DrawingClasses.so'))<0:
        raise RuntimeError('Rebuild the current plotting/DrawingClasses.so')
    ROOT.gInterpreter.AddIncludePath(str(REPO/'plotting'))
    if not ROOT.gInterpreter.Declare('#include "EECCard.h"\n#include "EECHistogramManager.h"\n#include "EECUnfoldConfiguration.h"'):
        raise RuntimeError('Cannot load existing EEC headers')
    if not hasattr(ROOT,'RooUnfoldBayes'): raise RuntimeError('RooUnfold dictionary unavailable')
    return ROOT


def error_enum(root):
    for owner in ('RooUnfolding','RooUnfold'):
        obj = getattr(root,owner,None)
        if obj is not None and hasattr(obj,'kCovariance'): return obj.kCovariance
    raise RuntimeError('Unsupported RooUnfold covariance enum')


def matrix_array(m):
    return np.array([[float(m[i][j]) for j in range(m.GetNcols())] for i in range(m.GetNrows())])


def t_matrix(root, values):
    values = covariance(values,check_psd=True)
    m = root.TMatrixD(*values.shape)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]): m[i][j] = float(values[i,j])
    return m


def unfolded_result(solver, mode):
    hm = next((n for n in ('Hunfold','Hreco') if hasattr(solver,n)),None)
    cm = next((n for n in ('Eunfold','Ereco') if hasattr(solver,n)),None)
    if hm is None or cm is None:
        raise RuntimeError('Unsupported result API; never fall back to GetMeasuredCov')
    h = getattr(solver,hm)(mode)
    if not h: raise RuntimeError('Solver returned no histogram')
    v = covariance(matrix_array(getattr(solver,cm)(mode)),h.GetNbinsX(),True)
    return h,v,dict(histogram_getter=hm,covariance_getter=cm)


def axis_edges(a):
    return np.array([a.GetBinLowEdge(i) for i in range(1,a.GetNbins()+2)])


def vector(h):
    return np.array([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)])


def variances(h):
    return np.array([h.GetBinError(i)**2 for i in range(1,h.GetNbinsX()+1)])


def hmatrix(h):
    return np.array([[h.GetBinContent(i,j) for j in range(1,h.GetNbinsY()+1)] for i in range(1,h.GetNbinsX()+1)])


def same_axis(a,b):
    a,b = axis_edges(a),axis_edges(b)
    return a.shape==b.shape and np.allclose(a,b,rtol=2e-6,atol=2e-8)


def clone(h, name):
    if not h: raise ValueError('Missing histogram: '+name)
    h = h.Clone(name); h.SetDirectory(0)
    return h


def check_flow(h):
    nx,ny = h.GetNbinsX(),h.GetNbinsY()
    points = [(0,),(nx+1,)] if h.GetDimension()==1 else [(i,j) for i in range(nx+2) for j in range(ny+2) if i in (0,nx+1) or j in (0,ny+1)]
    if any(h.GetBinContent(*p)!=0 or h.GetBinError(*p)!=0 for p in points):
        raise ValueError('Explicitly account for nonzero ROOT flow before this regular-guard-bin adapter: '+h.GetName())


def card_values(handle):
    d = handle.GetDirectory('JCard')
    if not d: raise ValueError('Missing JCard')
    result = {}
    for key in d.GetListOfKeys():
        name = key.GetName()
        if name in result: continue
        obj = d.Get(name)
        if obj.ClassName().startswith('TVectorT'):
            result[name] = [float(obj[i]) for i in range(obj.GetLwb(),obj.GetUpb()+1)]
        elif obj.InheritsFrom('TObjString'): result[name] = str(obj.GetString())
    return result


def equal(a,b):
    return np.array_equal(np.asarray(a,np.float32),np.asarray(b,np.float32)) if isinstance(a,list) and isinstance(b,list) else a==b


def assert_compatible(a,b,response_only=False):
    extra = ('McCorrelationType','MinJetPtUnfoldingReco','MinJetPtUnfoldingTruth','MinJetPtCut','MaxJetPtCut','SkipCovarianceMatrix','MatchJets') if response_only else ()
    for k in CONTRACT_KEYS+extra:
        if k not in a or k not in b: raise ValueError('Missing compatibility metadata: '+k)
        if not equal(a[k],b[k]): raise ValueError('Incompatible metadata: '+k)


class Product:
    def __init__(self, root, path, channel, tag, ordinary=False):
        self.path = Path(path).resolve(strict=True)
        self.handle = root.TFile.Open(str(self.path),'READ')
        if not self.handle or self.handle.IsZombie() or self.handle.TestBit(root.TFile.kRecovered):
            raise ValueError('Unreadable/recovered product: '+str(path))
        self.metadata = card_values(self.handle)
        fixed = dict(DataType=10 if channel['system']=='OO' else 12,JetRadius=.4,JetAxis=1,
            JetEtaCut=1.6,TrackEtaCut=2.4,MinTrackPtCut=2,MinPairTrackPt=2,SubeventCut=2,
            HFEnergyShift=0,JetPtWeight=0,SmearDeltaR=0,SmearEnergyWeight=0,BackgroundMethods=0)
        for k,v in fixed.items():
            if not equal(self.metadata.get(k),[v]): raise ValueError('Wrong pilot setting: '+k)
        if self.metadata.get('WeightExponent') != [float(channel['power'])]:
            raise ValueError('Use a separate explicit per-power projected product')
        self.card = root.EECCard(self.handle)
        self.ic = interval_index(self.metadata['CentralityBinEdges'],channel['activity'])
        self.it = physical_index(self.metadata['TrackPtBinEdgesEEC'],2.)
        if self.ic>=root.EECHistogramManager.kMaxCentralityBins or self.it>=self.card.GetNTrackPtBinsEEC():
            raise ValueError('Activity or threshold index exceeds native capacity')
        self.manager = root.EECHistogramManager(self.handle,self.card)
        self.manager.SetCentralityBinRange(self.ic,self.ic)
        self.manager.SetTrackPtBinRangeEEC(self.it,self.it)
        self.manager.SetJetPtBinRangeEEC(0,self.card.GetNJetPtBinsEEC()-1)
        self.manager.SetLoadJetPtUnfoldingHistograms(not ordinary)
        self.manager.SetLoadJetPtUnfoldingCovariance(True)
        self.manager.SetLoadEnergyEnergyCorrelators(ordinary)
        self.manager.LoadProcessedHistograms()
        self.jets = np.array(self.metadata['JetPtBinEdgesEEC'])
        self.tag = tag

    def response(self):
        m = clone(self.manager.GetHistogramJetPtUnfoldingMeasured(self.ic,self.it),self.tag+'_measured')
        t = clone(self.manager.GetHistogramJetPtUnfoldingTruth(self.ic,self.it),self.tag+'_truth')
        c = clone(self.manager.GetHistogramJetPtUnfoldingResponse(self.ic,self.it),self.tag+'_response')
        for h in (m,t,c): check_flow(h)
        if not same_axis(m.GetXaxis(),c.GetXaxis()) or not same_axis(t.GetXaxis(),c.GetYaxis()):
            raise ValueError('Response and marginal axes differ')
        return m,t,c

    def measured_cov(self,h):
        if self.metadata.get('McCorrelationType') != [1.]:
            raise ValueError('Generator-particle pilot requires covariance from RecoGen filling')
        c = self.manager.GetHistogramJetPtUnfoldingCovariance(0,self.ic,self.it)
        if not c: raise ValueError('Missing measured per-jet covariance, not replaced by Sumw2')
        if not same_axis(c.GetXaxis(),h.GetXaxis()) or not same_axis(c.GetYaxis(),h.GetXaxis()):
            raise ValueError('Measured covariance bin map mismatch')
        # Native LoadJetPtUnfoldingCovariance projects with normalizeToBinWidth=false.
        return covariance(hmatrix(c),h.GetNbinsX(),True)

    def ordinary(self,template,expected_mode):
        if self.metadata.get('McCorrelationType') != [float(expected_mode)]:
            raise ValueError('Wrong ordinary RecoGen/GenGen mode')
        n,dr = flattened_layout(axis_edges(template.GetXaxis()),self.jets)
        flat = clone(template,self.tag+'_bin_integrals'); flat.Reset()
        for j in range(len(self.jets)-1):
            h = self.manager.GetHistogramEnergyEnergyCorrelator(0,self.ic,j,self.it)
            if not h or not np.allclose(axis_edges(h.GetXaxis()),dr,rtol=2e-6,atol=2e-8):
                raise ValueError('Ordinary EEC missing or incompatible DeltaR bins')
            check_flow(h)
            for i in range(1,n+1):
                flat.SetBinContent(j*n+i,h.GetBinContent(i)*h.GetBinWidth(i))
                flat.SetBinError(j*n+i,h.GetBinError(i)*h.GetBinWidth(i))
        return flat

    def close(self):
        self.handle.Close()


def write_matrix(root,out,name,values,xedges,yedges):
    h = root.TH2D(name,name,len(xedges)-1,array('d',xedges),len(yedges)-1,array('d',yedges))
    h.SetDirectory(0)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]): h.SetBinContent(i+1,j+1,float(values[i,j]))
    out.cd(); h.Write(name)


def run_solver(root,response,measured,vm,iteration,name):
    solver = root.RooUnfoldBayes(response,measured,iteration)
    if hasattr(solver,'SetVerbose'): solver.SetVerbose(0)
    solver.SetMeasuredCov(t_matrix(root,vm))
    h,v,api = unfolded_result(solver,error_enum(root))
    return clone(h,name),v,api
