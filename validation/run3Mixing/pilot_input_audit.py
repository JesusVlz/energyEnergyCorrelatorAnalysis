#!/usr/bin/env python3
"""Bounded read-only schema/provenance and physical-identity audit for OO closure.

The event ledger, not a file-path split, defines the reproducible A/B partition.
Source ROOTs remain unchanged. Unknown embedding provenance is never certified.
"""
import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

SIGNAL_NAMESPACE = 'OO2025_embedded_existing_owner_pilot'
POOL_NAMESPACE = 'OO2025_HIJING_existing_owner_pool_pilot'
BACKGROUND_NAMESPACE = 'OO2025_HIJING_stored_kinematics_v1'
FILTERS = ('pprimaryVertexFilter', 'OOphfCoincFilterPF2Th4', 'pileupVertexFilter')
EVENT_BRANCHES = ('run', 'lumi', 'evt', 'vz', 'hiBin', 'hiHF_pf', 'weight')
GEN_BRANCHES = ('pt', 'eta', 'phi', 'chg', 'sube', 'pdg')
VETO_BRANCHES = ('nref', 'jtpt', 'rawpt', 'jteta', 'jtphi', 'WTAeta', 'WTAphi',
                 'jtPfCHF', 'jtPfNHF', 'jtPfCEF', 'jtPfNEF', 'jtPfMUF',
                 'jtPfCHM', 'jtPfNHM', 'jtPfCEM', 'jtPfNEM', 'jtPfMUM')
PRODUCER_ROOT = Path('/eos/user/j/jvelazqu/Foresting/CMSSW_15_0_11/src')
PRODUCER_FILES = ('HeavyIonsAnalysis/EventAnalysis/plugins/HiGenAnalyzer.cc',
                  'HeavyIonsAnalysis/EventAnalysis/python/HiGenAnalyzer_cfi.py',
                  'HeavyIonsAnalysis/Configuration/forest_miniAOD_run3_MC.py')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def physical_id(values):
    if len(values) != 3 or any(isinstance(x, bool) or not isinstance(x, int) or not 0 <= x < 1 << 64 for x in values):
        raise ValueError('run/lumi/event require three exact unsigned 64-bit integers')
    return tuple(values)


def event_key(namespace, values):
    """Byte-for-byte FNV identity used by src/OOEventIdentity.h."""
    physical_id(values)
    if not namespace:
        raise ValueError('Physical sample namespace required')
    value = 14695981039346656037
    for byte in namespace.encode()+b'\0'+b''.join(x.to_bytes(8, 'little') for x in values):
        value = ((value ^ byte)*1099511628211) & ((1 << 64)-1)
    return value


def background_key(signature):
    count, total, parity = signature.split(':')
    if int(count) <= 0:
        return None
    return event_key(BACKGROUND_NAMESPACE, [int(count), int(total,16), int(parity,16)])


def partition(values, seed, signature=None):
    if not 0 <= seed < 1 << 64:
        raise ValueError('Seed must fit unsigned 64-bit')
    key = background_key(signature) if signature else None
    return 'AB'[((event_key(POOL_NAMESPACE, values) if key is None else key) ^ seed) % 2]


def bounded_ranges(lengths, cap):
    """Global cap across all files; never a cap repeated for each file."""
    if cap <= 0:
        raise ValueError('Positive global scan cap required')
    remaining = cap
    result = []
    for entries in lengths:
        size = min(max(0, int(entries)), remaining)
        result.append(size)
        remaining -= size
    return result


def read_pool(path):
    files = [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
    if not files or len(set(files)) != len(files):
        raise ValueError('Pool list must be nonempty with no duplicate paths')
    return files


CPP = r'''
#include <vector>
#include <map>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
namespace oo_input_audit {
struct Summary {
 unsigned long long nh=0, particles=0, charged=0, backgroundParticles=0;
 std::uint64_t signatureSum=0,signatureXor=0;
 std::map<int,unsigned long long> labels;
 double minPt=std::numeric_limits<double>::infinity(),maxPt=0,maxAbsEta=0,maxAbsPhi=0;
};
inline std::uint64_t mix(std::uint64_t v){ v+=0x9e3779b97f4a7c15ULL;v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;v=(v^(v>>27))*0x94d049bb133111ebULL;return v^(v>>31); }
Summary summarize(const std::vector<float>& pt,const std::vector<float>& eta,const std::vector<float>& phi,
                  const std::vector<int>& chg,const std::vector<int>& sube,const std::vector<int>& pdg,bool signal){
 Summary out;
 const auto n=pt.size();
 if(eta.size()!=n||phi.size()!=n||chg.size()!=n||sube.size()!=n||pdg.size()!=n) throw std::runtime_error("Unaligned particle vectors");
 for(unsigned i=0;i<n;++i){
  if(!std::isfinite(pt[i])||!std::isfinite(eta[i])||!std::isfinite(phi[i])||pt[i]<0) throw std::runtime_error("Invalid particle coordinate");
  ++out.particles; ++out.labels[sube[i]]; out.charged+=(chg[i]!=0);
  out.minPt=std::min(out.minPt,double(pt[i]));out.maxPt=std::max(out.maxPt,double(pt[i]));
  out.maxAbsEta=std::max(out.maxAbsEta,double(std::abs(eta[i])));out.maxAbsPhi=std::max(out.maxAbsPhi,double(std::abs(phi[i])));
  const bool background=!signal || sube[i]>0;
  if(background){
   ++out.backgroundParticles;
   if(chg[i]!=0 && pt[i]>.7 && std::abs(eta[i])<2.4) ++out.nh;
   std::uint32_t a,b,c;std::memcpy(&a,&pt[i],4);std::memcpy(&b,&eta[i],4);std::memcpy(&c,&phi[i],4);
   auto h=mix(a)^mix(std::uint64_t(b)<<1)^mix(std::uint64_t(c)<<2)^mix(std::uint64_t(std::uint32_t(chg[i]))<<3);
   out.signatureSum+=h;out.signatureXor^=mix(h);
  }
 }
 return out;
}
}
'''


def branch_schema(tree):
    result = {}
    for branch in tree.GetListOfBranches():
        leaf = branch.GetListOfLeaves()[0] if len(branch.GetListOfLeaves()) else None
        result[branch.GetName()] = {'class': branch.GetClassName(),
                                    'type': leaf.GetTypeName() if leaf else None,
                                    'count_leaf': leaf.GetLeafCount().GetName() if leaf and leaf.GetLeafCount() else None}
    return result


def require_branches(tree, names):
    schema = branch_schema(tree)
    missing = sorted(set(names)-set(schema))
    if missing:
        raise ValueError(f'Missing branches in {tree.GetName()}: {missing}')
    tree.SetBranchStatus('*', 0)
    for name in names:
        tree.SetBranchStatus(name, 1)
    return schema


def inspect_file(ROOT, path):
    if not Path(path).is_file() or Path(path).stat().st_size <= 0:
        raise ValueError('Missing or empty input: '+path)
    file = ROOT.TFile.Open(path)
    if not file or file.IsZombie():
        raise ValueError('Cannot open input with ROOT: '+path)
    trees = {name: file.Get(name) for name in ('hiEvtAnalyzer/HiTree', 'HiGenParticleAna/hi', 'skimanalysis/HltTree', 'akCs4PFJetAnalyzer/t')}
    if not all(trees.values()):
        raise ValueError('Missing required native event/gen/skim/R4 tree: '+path)
    entries = int(trees['hiEvtAnalyzer/HiTree'].GetEntries())
    if any(int(t.GetEntries()) != entries for t in trees.values()):
        raise ValueError('Unaligned forest entries: '+path)
    schemas = {name: branch_schema(tree) for name, tree in trees.items()}
    for name, required in (('hiEvtAnalyzer/HiTree', EVENT_BRANCHES), ('HiGenParticleAna/hi', GEN_BRANCHES), ('skimanalysis/HltTree', FILTERS), ('akCs4PFJetAnalyzer/t', VETO_BRANCHES)):
        missing = set(required)-set(schemas[name])
        if missing:
            raise ValueError(f'Missing required branches {missing} in {path}:{name}')
    event_schema = schemas['hiEvtAnalyzer/HiTree']
    if any(event_schema[k]['type'] != 'UInt_t' for k in ('run', 'lumi')) or event_schema['evt']['type'] != 'ULong64_t':
        raise ValueError('Physical identity schema is not UInt_t/UInt_t/ULong64_t: '+path)
    info = file.Get('HiForestInfo/HiForest')
    info_text = {}
    if info and info.GetEntries():
        info.GetEntry(0)
        for branch in info.GetListOfBranches():
            name = branch.GetName()
            try:
                leaf = branch.GetListOfLeaves()[0]
                info_text[name] = ''.join(chr(int(leaf.GetValue(i)) & 255) for i in range(int(leaf.GetLen()))).split('\0',1)[0]
            except (AttributeError, TypeError):
                info_text[name] = 'UNREADABLE'
    metadata = {'path': path, 'size_bytes': Path(path).stat().st_size, 'root_uuid': str(file.GetUUID().AsString()),
                'source_entries': entries, 'schema': schemas, 'forest_info': info_text,
                'embedding_identity_fields': [name for name in schemas['hiEvtAnalyzer/HiTree'] if any(term in name.lower() for term in ('embed', 'background', 'hijing', 'parent'))],
                'reco_veto_required_branches_present': True,
                'final_state_convention': 'sta absent; local HiGenAnalyzer emits sta only when stableOnly=False; local cfi stableOnly=True' if 'sta' not in schemas['HiGenParticleAna/hi'] else 'sta branch available; must enforce status1'}
    file.Close()
    return metadata


def scan_file(ROOT, metadata, count, role, seed, writer, seen, signatures):
    file = ROOT.TFile.Open(metadata['path'])
    event, gen, skim = (file.Get(name) for name in ('hiEvtAnalyzer/HiTree', 'HiGenParticleAna/hi', 'skimanalysis/HltTree'))
    require_branches(event, EVENT_BRANCHES)
    gen_schema = require_branches(gen, GEN_BRANCHES)
    require_branches(skim, FILTERS)
    if 'sta' in gen_schema:
        raise ValueError('This bounded audit requires the source stableOnly convention; explicit status branch needs source-specific validation')
    namespace = SIGNAL_NAMESPACE if role == 'signal' else POOL_NAMESPACE
    counts, labels = Counter(), Counter()
    ranges = {'vz': [math.inf, -math.inf], 'hiBin': [math.inf, -math.inf], 'nh': [math.inf, -math.inf],
              'hf_raw': [math.inf, -math.inf], 'weight': [math.inf, -math.inf], 'particle_pt': [math.inf, -math.inf]}
    max_eta = max_phi = 0.
    for entry in range(count):
        event.GetEntry(entry); gen.GetEntry(entry); skim.GetEntry(entry)
        identity = physical_id([int(event.run), int(event.lumi), int(event.evt)])
        duplicate = identity in seen
        counts['duplicate_physical_ids'] += duplicate
        first = seen.setdefault(identity, [metadata['path'], entry])
        summary = ROOT.oo_input_audit.summarize(gen.pt, gen.eta, gen.phi, gen.chg, gen.sube, gen.pdg, role == 'signal')
        for label, value in summary.labels:
            labels[int(label)] += int(value)
        flags = {name: bool(getattr(skim, name)) for name in FILTERS}
        vz, hf, weight, hibin = float(event.vz), float(event.hiHF_pf), float(event.weight), int(event.hiBin)
        if any(not math.isfinite(v) for v in (vz, hf, weight)):
            raise ValueError('Nonfinite event coordinate/weight')
        pre_veto = abs(vz) <= 15 and all(flags.values())
        group = 'invalid' if not 0 <= hibin < 200 else ('90_100' if hibin >= 180 else '0_90')
        counts['scanned'] += 1; counts['source_hiBin_'+group] += 1
        counts['pre_veto_selected'] += pre_veto; counts['pre_veto_hiBin_'+group] += pre_veto
        counts['weight_negative'] += weight < 0; counts['weight_zero'] += weight == 0; counts['weight_unity'] += weight == 1
        counts['particles'] += int(summary.particles); counts['charged_particles'] += int(summary.charged)
        for name, passed in flags.items(): counts['passes_'+name] += passed
        nh = int(summary.nh)
        sig = f'{int(summary.backgroundParticles)}:{int(summary.signatureSum):016x}:{int(summary.signatureXor):016x}'
        if summary.backgroundParticles:
            signatures.setdefault(sig, []).append({'id': list(identity), 'source_file': metadata['path'], 'source_entry': entry})
        row = {'sample_namespace': namespace, 'id': list(identity), 'event_key': f'{event_key(namespace, identity):016x}',
               'source_file': metadata['path'], 'source_entry': entry,
               'vz': vz, 'hiBin': hibin, 'hf_raw': hf, 'nh': nh,
               'generator_weight': weight, 'filter_flags': flags, 'pre_veto_selection_pass': pre_veto,
               'reco_veto_pass': None, 'duplicate_physical_id': duplicate, 'first_physical_id_source': first,
               'background_identity': 'storedgen:'+sig if summary.backgroundParticles else None,
               'background_key': f'{background_key(sig):016x}' if summary.backgroundParticles else None,
               'background_identity_kind': 'OBSERVED_STORED_KINEMATICS_NOT_CERTIFIED_LINEAGE' if summary.backgroundParticles else 'UNKNOWN_EMPTY_BACKGROUND',
               'background_kinematics_signature': sig,
               'background_provenance': 'UNKNOWN_SHARED_EMBEDDING_ID',
               'pool': partition(identity, seed, sig) if role == 'pool' else None}
        writer(row)
        for name, value in (('vz',vz), ('hiBin',hibin), ('nh',nh), ('hf_raw',hf), ('weight',weight)):
            ranges[name][0] = min(ranges[name][0], value); ranges[name][1] = max(ranges[name][1], value)
        if summary.particles:
            ranges['particle_pt'][0] = min(ranges['particle_pt'][0], float(summary.minPt)); ranges['particle_pt'][1] = max(ranges['particle_pt'][1], float(summary.maxPt))
        max_eta = max(max_eta, float(summary.maxAbsEta)); max_phi = max(max_phi, float(summary.maxAbsPhi))
    file.Close()
    return {'counts': dict(counts), 'subevent_particle_counts': dict(labels), 'ranges': {k: [v if math.isfinite(v) else None for v in vs] for k, vs in ranges.items()},
            'max_abs_particle_eta': max_eta, 'max_abs_particle_phi': max_phi, 'entry_range_scanned': [0, count]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--signal', type=Path, required=True)
    parser.add_argument('--pool-list', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--signal-events', type=int, default=10000)
    parser.add_argument('--pool-events', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=20261030)
    parser.add_argument('--replace-derived', action='store_true', help='Regenerate the derived manifest/A/B lists; preserve old audit in its fresh result directory')
    args = parser.parse_args()
    if not 1 <= args.signal_events <= 10000 or not 1 <= args.pool_events <= 100000 or not 0 <= args.seed < 1 << 64:
        parser.error('Bounded audit requires signal 1..10000, global pool 1..100000, unsigned64 seed')
    if args.output.exists():
        parser.error('Fresh output directory required: '+str(args.output))
    derived = [args.pool_list.with_suffix('.manifest.json'), *[args.pool_list.with_name(args.pool_list.stem+'_'+pool+'.txt') for pool in 'AB']]
    if any(path.exists() for path in derived) and not args.replace_derived:
        parser.error('Derived manifest/A/B lists already exist; protect earlier audit outputs')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gInterpreter.Declare(CPP)
    start = time.monotonic()
    files = read_pool(args.pool_list)
    signal_meta = inspect_file(ROOT, str(args.signal))
    pool_meta = []
    for path in files:
        pool_meta.append(inspect_file(ROOT, path))
    args.output.mkdir(parents=True, exist_ok=False)
    counts = bounded_ranges([m['source_entries'] for m in pool_meta], args.pool_events)
    streams = {name: (args.output/(name+'.events.jsonl')).open('x') for name in ('signal', 'poolA', 'poolB')}
    seen = {'signal': {}, 'pool': {}}
    signatures = {'signal': {}, 'pool': {}}
    partition_counts = Counter()
    partition_files = {'A': set(), 'B': set()}
    def write(row):
        name = 'signal' if row['pool'] is None else 'pool'+row['pool']
        streams[name].write(json.dumps(row, allow_nan=False)+'\n')
        if row['pool'] is not None:
            partition_counts[row['pool']] += 1; partition_files[row['pool']].add(row['source_file'])
    try:
        signal_meta['audit'] = scan_file(ROOT, signal_meta, min(args.signal_events, signal_meta['source_entries']), 'signal', args.seed, write, seen['signal'], signatures['signal'])
        print('signal audit complete', flush=True)
        for metadata, count in zip(pool_meta, counts):
            if count:
                metadata['audit'] = scan_file(ROOT, metadata, count, 'pool', args.seed, write, seen['pool'], signatures['pool'])
                print('pool audit complete', metadata['path'], count, flush=True)
            else:
                metadata['audit'] = {'entry_range_scanned': [0,0], 'status': 'SCHEMA_ONLY_GLOBAL_CAP_REACHED'}
    finally:
        for stream in streams.values(): stream.close()
    overlap = sorted(set(seen['signal']) & set(seen['pool']))
    signature_overlap = sorted(set(signatures['signal']) & set(signatures['pool']))
    for pool, path in zip('AB', derived[1:]):
        with path.open('w' if args.replace_derived else 'x') as stream:
            stream.write('# Event-ledger partition required: this is not an independent whole-file pool.\n')
            stream.write('# Authoritative partition: '+str((args.output/('pool'+pool+'.events.jsonl')).resolve())+'\n')
            stream.write('\n'.join(path for path in files if path in partition_files[pool])+'\n')
    provenance_files = {str(PRODUCER_ROOT/path): digest(PRODUCER_ROOT/path) for path in PRODUCER_FILES if (PRODUCER_ROOT/path).exists()}
    manifest = {'status': 'BOUNDED_INPUT_AUDIT_NOT_PHYSICS_CLOSURE', 'seed': args.seed,
                'signal_global_event_cap': args.signal_events, 'pool_global_event_cap': args.pool_events,
                'signal': signal_meta, 'pool_files': pool_meta,
                'input_pool_list': str(args.pool_list.resolve()), 'input_pool_list_sha256': digest(args.pool_list),
                'command_argv': sys.argv,
                'source_commit': subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
                'root_version': ROOT.gROOT.GetVersion(), 'python_version': sys.version,
                'audit_tool_sha256': digest(__file__), 'identity_header_sha256': digest(Path(__file__).resolve().parents[2]/'src/OOEventIdentity.h'),
                'local_producer_files_sha256': provenance_files,
                'physical_id_source': 'hiEvtAnalyzer/HiTree: UInt_t run, UInt_t lumi, ULong64_t evt; no float conversion; generator/jet Int_t event fields unused',
                'physical_id_duplicates': {role: sum(m.get('audit',{}).get('counts',{}).get('duplicate_physical_ids',0) for m in ([signal_meta] if role=='signal' else pool_meta)) for role in ('signal','pool')},
                'signal_pool_raw_id_overlap_count': len(overlap), 'signal_pool_raw_id_overlap_examples': [list(identity) for identity in overlap[:20]],
                'background_signature_screen': {'method': 'order-independent count plus uint64 sum/xor of SplitMix64 combinations of stored float32 pt/eta/phi and charge; signal sube>0 versus all donor particles; diagnostic only',
                                                'exact_signature_overlap_count': len(signature_overlap), 'overlap_examples': signature_overlap[:20],
                                                'signal_events_with_nonempty_signature_overlap': sum(len(signatures['signal'][s]) for s in signature_overlap),
                                                'pool_events_with_nonempty_signature_overlap': sum(len(signatures['pool'][s]) for s in signature_overlap),
                                                'within_signal_repeated_nonempty_signatures': sum(len(v)-1 for v in signatures['signal'].values()),
                                                'within_pool_repeated_nonempty_signatures': sum(len(v)-1 for v in signatures['pool'].values()),
                                                'overlap_ledger': str((args.output/'background_identity_matches.jsonl').resolve()),
                                                'limitation': 'No matching signature does not exclude shared underlying HIJING; source configuration/campaign links unavailable'},
                'partition': {'rule': 'For nonempty background, (EventKey(BACKGROUND_NAMESPACE,count,signatureSum,signatureXor) XOR seed) modulo 2; for empty background EventKey(POOL_NAMESPACE,run,lumi,event) fallback. 0=A,1=B',
                              'background_namespace': BACKGROUND_NAMESPACE,
                              'namespace': POOL_NAMESPACE, 'events': dict(partition_counts), 'physical_id_disjoint': True, 'observed_nonempty_background_signature_disjoint': True,
                              'authoritative_event_ledgers': {name: {'path': str((args.output/(name+'.events.jsonl')).resolve()), 'sha256': digest(args.output/(name+'.events.jsonl'))} for name in streams},
                              'derived_file_lists': [str(path.resolve()) for path in derived[1:]],
                              'file_list_warning': 'A/B lists intentionally share source paths; selecting all file entries violates partition. Apply seeded observed-background signature fold (physical-ID fallback only for empty backgrounds) or authoritative ledger.'},
                'matching_definition': {'nh': 'stored stable charged generator pT>0.7 GeV and abs(eta)<2.4; signal sube>0 only; pure-HIJING all accepted particles by owner source provenance',
                                        'delta_vz_cm_max': .5, 'delta_nh_max': 'max(1,0.01*NHsignal)', 'independent_of_EEC_threshold': True,
                                        'OO_HF_transform': 'OFF; hf_raw copied unchanged; not a matching coordinate'},
                'provenance': {'pool_source': 'Owner-identified existing pure-HIJING HiForests; no new production needed',
                               'embedded_origin_labels': 'Local HiGenAnalyzer: sube0 packedGenParticlesSignal pointer membership; sube1 other packedGenParticles; -999 absent signal collection. Source labels audited in each scanned file.',
                               'stable_final_state': 'Supported by missing sta branch and local stableOnly=True producer; exact historical production config hash unavailable',
                               'source_particle_acceptance': 'Local HiGenAnalyzer cfi pT>=0.25 GeV, abs(eta)<=2.5; observed ranges recorded, not inferred untruncated spectrum',
                               'source_event_selection': 'OO pprimaryVertexFilter, OOphfCoincFilterPF2Th4, pileupVertexFilter, abs(vz)<=15 audited before Run3 reconstructed jet veto',
                               'reco_veto': 'Required full native R4 PF jet branches verified in every listed file; actual veto mask evaluated by analyzer, not this pre-veto audit',
                               'reco_veto_source_jet_retention': 'Current local producer config jetPtMin=15 GeV corroborates declared input retention. Exact historical file-production threshold not pinned by ROOT metadata; complete veto coverage remains conditional.',
                               'upstream_filter_loss_and_MB_representativeness': 'UNKNOWN; cannot infer from selected forest',
                               'shared_embedding_background_independence': 'UNKNOWN; no physical background provenance fields found. Exact raw physical IDs and stored kinematic signatures screened within bounded ranges only.',
                               'donor_weight_use': 'Raw stored weights audited; do not multiply donor hard-process weights into EEC terms'},
                'elapsed_seconds': time.monotonic()-start}
    with (args.output/'background_identity_matches.jsonl').open('x') as stream:
        for signature in signature_overlap:
            stream.write(json.dumps({'signature': signature, 'background_key': f'{background_key(signature):016x}', 'signals': signatures['signal'][signature], 'donors': signatures['pool'][signature]}, allow_nan=False)+'\n')
    for path in (args.output/'manifest.json', derived[0]):
        with path.open('w' if path == derived[0] and args.replace_derived else 'x') as stream: json.dump(manifest, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({'signal_scanned': signal_meta['audit']['counts']['scanned'], 'pool_scanned': sum(counts), 'pool_partition_counts': dict(partition_counts),
                      'raw_id_overlap': len(overlap), 'signature_overlap': len(signature_overlap), 'manifest': str(derived[0])}), flush=True)


if __name__ == '__main__':
    main()
