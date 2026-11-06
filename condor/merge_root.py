#!/usr/bin/env python3
"""Merge EEC ROOT files, retaining configuration and checking every histogram bin.

The campaign worker owns artifact hashes, ledger checks, and publication. This
module requires PyROOT and hadd and checks the ROOT files in local scratch.
"""
import json
import math
from pathlib import Path
import subprocess
import tempfile

RUNTIME_CARD_KEYS = frozenset(('InputFiles', 'MixingFileList', 'MixingFileListPath',
                               'OOClosureOutputPrefix'))
_CPP_READY = False


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _root():
    import ROOT
    global _CPP_READY
    ROOT.gROOT.SetBatch(True)
    if not _CPP_READY:
        ROOT.gInterpreter.Declare(r'''
#include <TH1.h>
#include <THnBase.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
namespace EECMergeRoot {
struct Bin { long double value=0, absolute=0, error2=0; };
struct Histogram { std::unordered_map<std::string,Bin> bins; long double entries=0; bool errors=false; };
inline void finite(double x) { if (!std::isfinite(x)) throw std::runtime_error("Nonfinite histogram value"); }
std::vector<double> Summary(TObject* object) {
  auto* one=dynamic_cast<TH1*>(object); auto* sparse=dynamic_cast<THnBase*>(object);
  if (!one && !sparse) throw std::runtime_error("Unsupported histogram");
  const double entries=one ? one->GetEntries() : sparse->GetEntries(); finite(entries);
  if (entries<0) throw std::runtime_error("Negative histogram entries");
  const auto n=one ? static_cast<Long64_t>(one->GetNcells()) : sparse->GetNbins();
  long double sum=0,absolute=0,error2=0;
  for (Long64_t i=0;i<n;++i) {
    const double v=one ? one->GetBinContent(i) : sparse->GetBinContent(i);
    const double e=one ? std::pow(one->GetBinError(i),2) : sparse->GetBinError2(i);
    finite(v); finite(e); if (e<0) throw std::runtime_error("Negative histogram error square");
    sum+=v; absolute+=std::abs(v); error2+=e;
  }
  return {entries,double(sum),double(absolute),double(error2),double(n)};
}
class Accumulator {
  std::map<std::string,Histogram> histograms;
public:
  void Add(const std::string& path,TObject* object) {
    auto* one=dynamic_cast<TH1*>(object); auto* sparse=dynamic_cast<THnBase*>(object);
    if (!one && !sparse) throw std::runtime_error("Unsupported histogram "+path);
    auto& h=histograms[path];
    h.entries+=one ? one->GetEntries() : sparse->GetEntries();
    h.errors=one ? one->GetSumw2N()>0 : sparse->GetCalculateErrors();
    const auto n=one ? static_cast<Long64_t>(one->GetNcells()) : sparse->GetNbins();
    std::vector<Int_t> coordinate(sparse ? sparse->GetNdimensions() : 1);
    for (Long64_t i=0;i<n;++i) {
      const double value=one ? one->GetBinContent(i) : sparse->GetBinContent(i,coordinate.data());
      const double e=one ? std::pow(one->GetBinError(i),2) : sparse->GetBinError2(i);
      finite(value); finite(e);
      const std::string key=one ? std::to_string(i) : std::string(reinterpret_cast<const char*>(coordinate.data()),coordinate.size()*sizeof(Int_t));
      auto& b=h.bins[key]; b.value+=value; b.absolute+=std::abs(value); b.error2+=e;
    }
  }
  void Check(const std::string& path,TObject* object) const {
    const auto found=histograms.find(path);
    if (found==histograms.end()) throw std::runtime_error("Unexpected histogram "+path);
    const auto& h=found->second;
    auto* one=dynamic_cast<TH1*>(object); auto* sparse=dynamic_cast<THnBase*>(object);
    const double entries=one ? one->GetEntries() : sparse->GetEntries();
    if (!std::isfinite(entries) || std::abs(entries-h.entries)>std::max(1.e-6L,std::abs(h.entries)*1.e-12L))
      throw std::runtime_error("Histogram entries do not sum: "+path);
    const auto n=one ? static_cast<Long64_t>(one->GetNcells()) : sparse->GetNbins();
    std::vector<Int_t> coordinate(sparse ? sparse->GetNdimensions() : 1);
    std::unordered_map<std::string,bool> seen;
    for (Long64_t i=0;i<n;++i) {
      const double v=one ? one->GetBinContent(i) : sparse->GetBinContent(i,coordinate.data());
      const double e=one ? std::pow(one->GetBinError(i),2) : sparse->GetBinError2(i);
      finite(v); finite(e);
      const std::string key=one ? std::to_string(i) : std::string(reinterpret_cast<const char*>(coordinate.data()),coordinate.size()*sizeof(Int_t));
      const auto expected=h.bins.find(key);
      if (expected==h.bins.end()) {
        if (v!=0 || (h.errors && e!=0)) throw std::runtime_error("Unexpected populated bin: "+path);
        continue;
      }
      const auto& b=expected->second; seen.emplace(key,true);
      // Float sums round at each addition. Magnitudes also bound cancellation.
      const long double tolerance=path=="nEvents" ? 1.e-6L : std::max(1.e-7L,b.absolute*5.e-5L);
      if (std::abs(v-b.value)>tolerance) throw std::runtime_error("Histogram bin does not sum: "+path);
      if (h.errors && std::abs(e-b.error2)>std::max(1.e-7L,std::abs(b.error2)*5.e-5L))
        throw std::runtime_error("Histogram squared errors do not sum: "+path);
    }
    for (const auto& item:h.bins) if (!seen.count(item.first) &&
        (item.second.value!=0 || (h.errors && item.second.error2!=0)))
      throw std::runtime_error("Missing populated merged bin: "+path);
  }
};
}
''')
        _CPP_READY = True
    return ROOT


def _open(path, mode='READ'):
    ROOT = _root()
    file = ROOT.TFile.Open(str(path), mode)
    require(bool(file) and file.IsOpen() and not file.IsZombie() and
            not file.TestBit(ROOT.TFile.kRecovered), 'Unreadable/zombie/recovered ROOT: '+str(path))
    return file


def _unique_keys(directory):
    keys = {}
    for key in directory.GetListOfKeys():
        name = str(key.GetName())
        require(name not in keys, 'Multiple ROOT key cycles: '+directory.GetPath()+'/'+name)
        keys[name] = key
    return keys


def _read_key(directory, key):
    object = key.ReadObj()
    # THnSparse/vector/TNamed ReadObj results are not registered with the
    # directory and otherwise leak across hundreds of inputs in PyROOT.
    if not directory.GetList().FindObject(object):
        _root().SetOwnership(object, True)
    return object


def _card_value(object):
    if object.InheritsFrom('TObjString'):
        return {'class': str(object.ClassName()), 'value': str(object.GetString())}
    if object.InheritsFrom('TVectorT<float>') or object.InheritsFrom('TVectorT<double>'):
        # operator() preserves the one-based lower bound of ConfigurationCard.
        values = [float(object(i)) for i in range(object.GetLwb(), object.GetUpb()+1)]
        require(all(math.isfinite(value) for value in values), 'Nonfinite numeric card entry')
        return {'class': str(object.ClassName()), 'lower_bound': int(object.GetLwb()), 'value': values}
    raise ValueError('Unsupported JCard object: '+str(object.ClassName()))


def _axis(axis):
    return {'edges': [float(axis.GetBinLowEdge(i)) for i in range(1, axis.GetNbins()+2)],
            'labels': [str(axis.GetBinLabel(i)) for i in range(1, axis.GetNbins()+1)]}


def _walk(directory, prefix=''):
    for name, key in _unique_keys(directory).items():
        if not prefix and name in ('JCard', 'MergeProvenance'):
            continue
        object = _read_key(directory, key)
        path = prefix+name
        if object.InheritsFrom('TDirectory'):
            yield from _walk(object, path+'/')
        else:
            yield path, object


def inspect_root(path, _accumulator=None):
    """Return JSON-safe card, histogram schema/totals and exposure metadata.

    Check every dense or allocated sparse bin including under/overflow for
    finite contents and errors. Unsupported statistical objects fail closed.
    """
    ROOT = _root()
    file = _open(path)
    try:
        directory = file.GetDirectory('JCard')
        require(bool(directory), 'Missing JCard: '+str(path))
        card = {name: _card_value(_read_key(directory,key)) for name, key in _unique_keys(directory).items()}
        histograms, metadata, runtime_metadata = {}, {}, {}
        for name, object in _walk(file):
            is_one, is_sparse = object.InheritsFrom('TH1'), object.InheritsFrom('THnBase')
            if is_one or is_sparse:
                require(not object.InheritsFrom('TProfile') and not object.InheritsFrom('TProfile2D') and
                        not object.InheritsFrom('TProfile3D'), 'Profile merging is unsupported: '+name)
                summary = list(ROOT.EECMergeRoot.Summary(object))
                axes = ([object.GetAxis(i) for i in range(object.GetNdimensions())] if is_sparse else
                        [object.GetXaxis(), object.GetYaxis(), object.GetZaxis()][:object.GetDimension()])
                histograms[name] = {'class': str(object.ClassName()), 'axes': [_axis(axis) for axis in axes],
                    'errors_stored': bool(object.GetCalculateErrors() if is_sparse else object.GetSumw2N()),
                    **dict(zip(('entries','sum_content','sum_absolute_content','sum_error2','allocated_bins'), summary))}
                if _accumulator is not None:
                    _accumulator.Add(name, object)
            else:
                require(str(object.ClassName()) in ('TNamed','TObjString'),
                        'Unsupported ROOT object '+name+': '+str(object.ClassName()))
                value = {'class': str(object.ClassName()), 'value': str(object.GetString()) if
                         object.InheritsFrom('TObjString') else str(object.GetTitle())}
                (runtime_metadata if name == 'InputFiles' else metadata)[name] = value
        for required in ('nEvents', 'inclusiveJet', 'energyEnergyCorrelator'):
            require(required in histograms, 'Missing histogram '+required+': '+str(path))
        events = file.Get('nEvents')
        require(events.InheritsFrom('TH1') and events.GetDimension() == 1, 'nEvents is not one-dimensional')
        require(events.GetNbinsX() >= 9, 'nEvents is missing event-cut bins')
        return {'path': str(path), 'card': card,
                'invariant_card': {key:value for key,value in card.items() if key not in RUNTIME_CARD_KEYS},
                'ignored_card': {key:value for key,value in card.items() if key in RUNTIME_CARD_KEYS},
                'metadata': metadata, 'runtime_metadata': runtime_metadata, 'histograms': histograms,
                'all_events': float(events.GetBinContent(1)),
                'inclusive_jet_entries': histograms['inclusiveJet']['entries'],
                'eec_entries': histograms['energyEnergyCorrelator']['entries'],
                'n_events_bins': [float(events.GetBinContent(i)) for i in range(events.GetNcells())]}
    finally:
        file.Close()


def _schema(summary):
    return {name: {key: value[key] for key in ('class', 'axes', 'errors_stored')}
            for name,value in summary['histograms'].items()}


def _restore_card(first_input, output, provenance):
    ROOT = _root()
    source, destination = _open(first_input), _open(output, 'UPDATE')
    try:
        destination.Delete('JCard;*')
        target = destination.mkdir('JCard')
        require(bool(target), 'Could not recreate JCard')
        target.cd()
        for name,key in _unique_keys(source.GetDirectory('JCard')).items():
            require(_read_key(source.GetDirectory('JCard'),key).Write(name, ROOT.TObject.kOverwrite)>0,
                    'Could not restore JCard/'+name)
        # hadd creates cycles for non-additive metadata too. Equality was
        # checked across inputs, so retain one canonical metadata object.
        for path, object in _walk(source):
            if object.InheritsFrom('TH1') or object.InheritsFrom('THnBase'):
                continue
            parent, _, name = path.rpartition('/')
            target = destination.GetDirectory(parent) if parent else destination
            require(bool(target), 'Missing metadata directory '+parent)
            target.Delete(name+';*')
            target.cd()
            require(object.Write(name,ROOT.TObject.kOverwrite)>0, 'Could not restore metadata '+path)
        destination.cd()
        # Each merged parent carries its own provenance. Remove all cycles
        # created by hadd before writing this reduction's canonical record.
        destination.Delete('MergeProvenance;*')
        require(ROOT.TObjString(json.dumps(provenance,sort_keys=True,allow_nan=False)).Write(
            'MergeProvenance',ROOT.TObject.kOverwrite)>0, 'Could not write merge provenance')
    finally:
        destination.Close()
        source.Close()


def merge_root(inputs, output, expected_totals=None):
    """hadd without skipping, restore the first card, validate every merged bin.

    Inputs must be unique, output must not exist. Optional expected_totals keys:
    all_events, inclusive_jet_entries, eec_entries, n_events_bins, from verified
    production records. Per-input runtime card paths remain in MergeProvenance.
    """
    ROOT = _root()
    paths = [str(Path(path).resolve()) for path in inputs]
    output = Path(output).resolve()
    require(paths and len(paths)==len(set(paths)), 'Empty or duplicate merge input paths')
    require(str(output) not in paths and not output.exists(), 'Merge output already exists or aliases an input')
    require(all('\n' not in path and '\r' not in path for path in paths), 'Newline in input path')
    accumulator = ROOT.EECMergeRoot.Accumulator()
    records = []
    for path in paths:
        record = inspect_root(path, _accumulator=accumulator)
        if records:
            require(record['invariant_card']==records[0]['invariant_card'], 'Invariant JCard differs: '+path)
            require(_schema(record)==_schema(records[0]), 'Histogram schema differs: '+path)
            require(record['metadata']==records[0]['metadata'], 'ROOT metadata differs: '+path)
        records.append(record)
    totals = {key: math.fsum(record[key] for record in records)
              for key in ('all_events','inclusive_jet_entries','eec_entries')}
    totals['n_events_bins'] = [math.fsum(values) for values in zip(*(record['n_events_bins'] for record in records))]
    for key,value in (expected_totals or {}).items():
        require(key in totals, 'Unsupported expected total: '+key)
        require(value==totals[key], 'ROOT inputs differ from verified production '+key)
    provenance = {'input_count':len(paths), 'canonical_card_source':paths[0],
        'runtime_card_keys_excluded_from_equality':sorted(RUNTIME_CARD_KEYS),
        'inputs':[{'path':record['path'], 'runtime_card':record['ignored_card'],
                   'runtime_metadata':record['runtime_metadata']} for record in records]}
    with tempfile.NamedTemporaryFile(mode='w',prefix='merge-root-',suffix='.inputs.txt',dir=output.parent) as source_list:
        source_list.write('\n'.join(paths)+'\n'); source_list.flush()
        command = ['hadd','-n','20',str(output),'@'+source_list.name]
        subprocess.run(command, check=True)
    _restore_card(paths[0],output,provenance)
    merged = inspect_root(output)
    require(merged['card']==records[0]['card'], 'Canonical merged JCard differs from first input')
    require(_schema(merged)==_schema(records[0]), 'Merged histogram schema differs')
    require(merged['metadata']==records[0]['metadata'], 'Merged ROOT metadata differs')
    file = _open(output)
    try:
        for name,object in _walk(file):
            if name in merged['histograms']:
                accumulator.Check(name,object)
    finally:
        file.Close()
    for key,value in totals.items():
        require(merged[key]==value, 'Merged ROOT exposure does not sum: '+key)
    return {'status':'PASS_ALL_HISTOGRAM_BINS_AND_CANONICAL_CARD', 'input_count':len(paths),
            'command':command, 'expected_totals':totals, 'root_validation':merged,
            'provenance':provenance,
            'validation_scope':'All histogram entries and dense/allocated sparse bin contents; stored squared errors; axes, labels, classes, invariant card and metadata; exact nEvents bins and exposures.',
            'float_bin_relative_tolerance':5e-5}
