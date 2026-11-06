#!/usr/bin/env python3
"""Strict pilot metadata gates; never perform unfolding or claim physics closure."""
import math
import struct

def float32(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Nonfinite physical value')
    return struct.unpack('f', struct.pack('f', value))[0]

def edge_indices(requested, available):
    """Exact equality after JCard Float_t serialization, never nearest-bin mapping."""
    stored = [float32(x) for x in available]
    if len(stored) < 2 or any(a >= b for a,b in zip(stored, stored[1:])):
        raise ValueError('Invalid or serialization-degenerate available edges')
    wanted = [float32(x) for x in requested]
    if any(a >= b for a,b in zip(wanted, wanted[1:])):
        raise ValueError('Requested edges are not increasing')
    try:
        return [stored.index(x) for x in wanted]
    except ValueError as error:
        raise ValueError('Unsupported physical edge; refill raw forests') from error

def read_card(path):
    card = {}
    for line in path.read_text().splitlines():
        fields = line.split('#', 1)[0].split()
        if not fields: continue
        if fields[0] in card: raise ValueError('Duplicate key: '+fields[0])
        try: card[fields[0]] = [float(x) for x in fields[1:]]
        except ValueError: card[fields[0]] = fields[1:]
    return card

SELECTIONS = dict(JetType=1, JetAxis=1, JetRadius=.4, ChannelRadiusCode=4,
    McCorrelationType=1, MatchJets=0, TriggerSelection=0,
    MinTrackPtCut=2, MinPairTrackPt=2, TrackEtaCut=2.4, MaxTrackPtCut=500,
    JetEtaCut=1.6, MinJetPtCut=40, MaxJetPtCut=5020,
    MinJetPtUnfoldingReco=40, MinJetPtUnfoldingTruth=10,
    Run3JetID=1, ApplyJetVeto=1, JetUncertainty=0, JetPtWeight=0,
    SmearDeltaR=0, SmearEnergyWeight=0, DisableTrackPairEfficiencyCorrection=1,
    BackgroundMethods=0, HFEnergyShift=0, SkipCovarianceMatrix=0,
    FilledHistograms=211, SubeventCut=2, MaxEventsPerFile=10000)
EDGES = [10,20,40,80,120,160,200,300,500,5020]

def validate_card(card, system, power, activity_edges):
    expected = dict(SELECTIONS, DataType=10 if system=='OO' else 12)
    for key,value in expected.items():
        if key not in card or len(card[key])!=1 or float32(card[key][0])!=float32(value):
            raise ValueError('Incompatible pilot selection: '+key)
    if card.get('WeightExponent') != [float(power)] or power not in (1,2):
        raise ValueError('Absent or incompatible energy power')
    if [float32(x) for x in card.get('CentralityBinEdges',[])] != [float32(x) for x in activity_edges]:
        raise ValueError('Incompatible activity selection')
    edges=card.get('JetPtBinEdgesEEC',[])
    if len(edges)!=len(EDGES) or edge_indices(EDGES,edges)!=list(range(len(EDGES))):
        raise ValueError('Unexpected pilot support')
    first,last=edge_indices([40,160],edges)
    if card.get('FirstUnfoldedJetPtBin') != [first] or card.get('LastUnfoldedJetPtBin') != [last-1]:
        raise ValueError('Output selector is not exactly 40--160 GeV')
    edge_indices([2],card.get('TrackPtBinEdgesEEC',[]))
    return True

def validate_product_metadata(expected, actual, power):
    """Extra cut gate supplementing the unchanged EECUnfoldConfiguration interface."""
    keys=set(SELECTIONS)-{'MaxEventsPerFile'}
    keys.update(('DataType','FirstUnfoldedJetPtBin','LastUnfoldedJetPtBin',
        'JetVetoInputMinPt','ZVertexCut','LowPtHatCut','HighPtHatCut','MultiplicityMode'))
    for key in sorted(keys):
        if key not in actual or [float32(x) for x in actual[key]] != [float32(x) for x in expected[key]]:
            raise ValueError('Missing/incompatible embedded selection: '+key)
    if float(power) not in actual.get('WeightExponent',[]):
        raise ValueError('Requested power absent from embedded card')
    for key in ('CentralityBinEdges','JetPtBinEdgesEEC','TrackPtBinEdgesEEC'):
        indices=edge_indices(expected[key],actual.get(key,[]))
        if indices != list(range(len(actual[key]))):
            raise ValueError('Incompatible full support: '+key)
    return True
