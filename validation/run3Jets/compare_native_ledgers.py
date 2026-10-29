#!/usr/bin/env python3
"""Unity-weight R4 payload hypotheses on retained native branch-audit ledgers.

No forest is opened. These kinematic-only samples do not reproduce analyzer
selection, event weights, PF Jet ID, vetoes or EEC closure.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from compare_calibrations import Corrector, PT_EDGES, REPO, native_match, setup_library, weighted_profile

SAMPLES = {'OO_data_production': ('OO', True), 'OO_MC_local': ('OO', False),
           'ppRef_data_production': ('pp', True), 'ppRef_MC_local': ('pp', False)}
RECO_FLOATS = ('rawpt', 'jtpt', 'jteta', 'jtphi', 'WTAeta', 'WTAphi', 'jtarea',
               'jtPfCHF', 'jtPfNHF', 'jtPfCEF', 'jtPfNEF', 'jtPfMUF')
RECO_INTS = ('jtPfCHM', 'jtPfNHM', 'jtPfCEM', 'jtPfNEM', 'jtPfMUM')
REF_FLOATS = ('refpt', 'refeta', 'refphi')
GEN_FLOATS = ('genpt', 'geneta', 'genphi')
LIMITS = {
    'selection': 'rawpt > 0 and abs(WTAeta) < 1.6 only; original reco indices retained',
    'weights': 'unity diagnostic weights; generator/vz/activity/cross-section weights NOT_RUN',
    'missing_selection_inputs': 'vz, hiBin and trackMax absent; full event/jet quality and activity splits NOT_RUN',
    'data_selection': 'trigger, PF Jet ID, veto-map applicability and event veto NOT_RUN',
    'candidate': 'NOT_APPROVED raw-input eta/pT-only L2 hypothesis; pp V2 CHS/producer and includes-L1 prescription unresolved; no L1/rho/area invented',
    'area': 'stored jtarea retained verbatim, including OO zeros; no physical-area validation or substitution',
    'matching': 'legacy native first ref/gen match; unmatched reco jets retained in pT comparisons; response only for positive matched gen pT',
    'physics': 'bounded jet-level diagnostic only; EEC n=1/n=2 closure, JER, physics closure and production approval NOT_RUN',
    'coverage': 'production PSet and veto-jet coverage above 15 GeV unresolved; pp data reference threshold 20 GeV',
    'radius': 'native R4 only; R8 remains independent raw-only work',
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact_integer(value, name, maximum=None):
    if type(value) is not int or value < 0 or (maximum is not None and value > maximum):
        raise ValueError(name+' must be an exact nonnegative JSON integer in range')
    return value


def validate_event(row, data):
    """Validate every native array before any selection, preserving its ordering."""
    if not isinstance(row, dict) or row.get('original_order_preserved') is not True:
        raise ValueError('original_order_preserved must be true')
    exact_integer(row.get('entry'), 'entry')
    v = row.get('values')
    if not isinstance(v, dict):
        raise ValueError('values must be an object')
    nref = exact_integer(v.get('nref'), 'nref')
    for key, maximum in [('run', 2**32-1), ('lumi', 2**32-1), ('evt', 2**64-1)]:
        exact_integer(v.get(key), key, maximum)
    arrays = [(key, nref, False) for key in RECO_FLOATS]+[(key, nref, True) for key in RECO_INTS]
    if not data:
        ngen = exact_integer(v.get('ngen'), 'ngen')
        arrays += [(key, nref, False) for key in REF_FLOATS]+[(key, ngen, False) for key in GEN_FLOATS]
    for key, count, integer in arrays:
        a = v.get(key)
        if not isinstance(a, list) or len(a) != count:
            raise ValueError(key+' array length must equal '+('ngen' if key in GEN_FLOATS else 'nref'))
        for value in a:
            if integer:
                exact_integer(value, key)
            elif type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError(key+' must contain finite numeric values')
    return v


def verify_artifact(evidence, manifest, name):
    p = evidence/name
    record = manifest.get(name)
    if not isinstance(record, dict) or record.get('sha256') != digest(p):
        raise ValueError('Evidence manifest hash mismatch or missing entry: '+name)
    if record.get('bytes') != p.stat().st_size:
        raise ValueError('Evidence manifest byte count mismatch: '+name)
    return {'path': str(p.resolve()), **record}


def load_sample(evidence, manifest, name, data):
    provenance = {suffix: verify_artifact(evidence, manifest, name+suffix) for suffix in ('.json', '.events.jsonl')}
    meta = json.loads((evidence/(name+'.json')).read_text())
    tree = 'akCs4PFJetAnalyzer/t' if SAMPLES[name][0] == 'OO' else 'ak4PFJetAnalyzer/t'
    if meta.get('sample') != name or meta.get('status') != 'INSPECTED':
        raise ValueError('Sample metadata identity/status mismatch: '+name)
    if meta.get('expected_tree') != tree or meta.get('event_tree') != 'hiEvtAnalyzer/HiTree':
        raise ValueError('Native R4/event tree metadata mismatch: '+name)
    if meta.get('array_length_errors') != []:
        raise ValueError('Retained input reports array-length errors: '+name)
    spec = [('nref', 'Int_t', None), ('run', 'UInt_t', None), ('lumi', 'UInt_t', None), ('evt', 'ULong64_t', None)]
    spec += [(k, 'Float_t', 'nref') for k in RECO_FLOATS]+[(k, 'Int_t', 'nref') for k in RECO_INTS]
    if not data:
        spec += [('ngen', 'Int_t', None)]+[(k, 'Float_t', 'nref') for k in REF_FLOATS]+[(k, 'Float_t', 'ngen') for k in GEN_FLOATS]
    for key, dtype, count in spec:
        b = meta.get('branches', {}).get(key, {})
        expected_tree = 'hiEvtAnalyzer/HiTree' if key in ('run', 'lumi', 'evt') else tree
        if not b.get('exists') or b.get('type') != dtype or b.get('count_branch') != count or b.get('tree') != expected_tree:
            raise ValueError('Native branch type/count/tree mismatch: '+key)
    rows, identities = [], set()
    for line in (evidence/(name+'.events.jsonl')).read_text().splitlines():
        row = json.loads(line)
        v = validate_event(row, data)
        if row['entry'] != len(rows):
            raise ValueError('Source entries must be contiguous in retained native order')
        identity = tuple(v[k] for k in ('run', 'lumi', 'evt'))
        if identity in identities:
            raise ValueError('Duplicate native event identity')
        identities.add(identity)
        rows.append(row)
    if len(rows) != meta.get('events_inspected'):
        raise ValueError('Retained event count differs from metadata')
    if sum(r['values']['nref'] for r in rows) != meta.get('jet_values_inspected_nref_sum'):
        raise ValueError('Retained reco count differs from metadata')
    return rows, meta, provenance


def payload_chains(system, data, upstream, repo=REPO):
    """Only pp data appends residuals, with separate L2-only/residual-only changes."""
    folder = repo/'jetEnergyCorrections'
    old = folder/('L2Relative_AK4PF_OO_V1.txt' if system == 'OO' else 'Prompt24HIpp_V1_MC_L2Relative_AK4PF.txt')
    new = upstream/('txt/2025OO/Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt' if system == 'OO' else 'txt/2024ppRef/Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt')
    chains = {'installed': [old], 'candidate_raw_input': [new]}
    if system == 'pp' and data:
        a = folder/'Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt'
        b = upstream/'txt/2024ppRef/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt'
        chains = {'installed': [old, a], 'candidate_raw_input': [new, b],
                  'residual_only_candidate': [old, b], 'l2_only_candidate': [new, a]}
    return chains


def compare_events(events, system, data, correctors):
    rows = []
    counts = Counter(events=len(events), native_reco=sum(e['values']['nref'] for e in events),
                     selected_reco=0, selected_stored_zero_area=0)
    if not data:
        counts.update(matched_reco=0, unmatched_reco=0, response_reco=0)
    for event in events:
        v = event['values']
        jets = SimpleNamespace(**v)
        for index in range(v['nref']):
            raw, eta, phi = (v[k][index] for k in ('rawpt', 'WTAeta', 'WTAphi'))
            if raw <= 0 or abs(eta) >= 1.6:
                continue
            match = None if data else native_match(jets, index)
            genpt = None if match is None else v['genpt'][match]
            row = {'id': [v[k] for k in ('run', 'lumi', 'evt')], 'source_entry': event['entry'],
                   'jet_index': index, 'system': system, 'data': data, 'diagnostic_weight': 1,
                   'gen_index': match, 'genpt': genpt, 'geneta': None if match is None else v['geneta'][match],
                   'rawpt': raw, 'stored_jtpt': v['jtpt'][index], 'WTAeta': eta, 'WTAphi': phi,
                   'jteta': v['jteta'][index], 'jtphi': v['jtphi'][index], 'stored_jtarea': v['jtarea'][index], 'response': {}}
            for name, corrector in correctors.items():
                value, clamped = corrector.evaluate(raw, eta, phi)
                row[name] = value
                row['response'][name] = value/genpt if genpt is not None and genpt > 0 else None
                counts[name+'_first_level_clamped'] += int(clamped)
            row['candidate_over_installed'] = row['candidate_raw_input']/row['installed']
            rows.append(row)
            counts['selected_reco'] += 1
            if not data:
                counts['matched_reco' if match is not None else 'unmatched_reco'] += 1
                counts['response_reco'] += int(genpt is not None and genpt > 0)
            counts['selected_stored_zero_area'] += int(row['stored_jtarea'] == 0)
    return rows, dict(counts)


def save_plot(fig, stem):
    fig.tight_layout()
    paths = []
    for suffix in ('.png', '.pdf'):
        path = stem.with_suffix(suffix)
        fig.savefig(path, dpi=150)
        paths.append(str(path))
    plt.close(fig)
    return paths


def profile(rows, coordinate, edges, numerator, denominator):
    result = weighted_profile(np.array([r[coordinate] for r in rows]),
        np.array([r[numerator]/r[denominator] for r in rows]), np.ones(len(rows)), edges)
    for cell in result:
        cell['sparse'] = cell['count'] < 10
    return result


def draw_profile(axis, cells, statistic, label):
    axis.plot([(p['low']+p['high'])/2 for p in cells],
              [p[statistic] if p[statistic] is not None else np.nan for p in cells], '.-', label=label)


def plot_sample(rows, name, data, chains, output):
    variants, profiles = list(chains), {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    fig.suptitle(name+': unity weights, kinematic-only; candidates unapproved', fontsize=10)
    for key in ['rawpt', 'stored_jtpt', *variants]:
        axes[0].hist([r[key] for r in rows], bins=PT_EDGES, histtype='step', label=key)
    axes[0].set(xlabel='jet pT [GeV]', ylabel='retained jets')
    axes[0].legend(fontsize=7)
    for key in variants[1:]:
        axes[1].hist([r[key]/r['installed'] for r in rows], bins=np.linspace(.5, 1.5, 61), histtype='step', label=key)
    axes[1].set(xlabel='variant / installed', ylabel='retained jets')
    axes[1].legend(fontsize=7)
    plots = save_plot(fig, output/(name+'_spectra_and_ratios'))
    if data and 'residual_only_candidate' in variants:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        fig.suptitle(name+': isolated residual/L2 changes; candidates unapproved', fontsize=10)
        for key in variants[1:]:
            for axis, coordinate, edges in [(axes[0], 'installed', PT_EDGES), (axes[1], 'WTAeta', np.linspace(-1.6, 1.6, 9))]:
                cells = profile(rows, coordinate, edges, key, 'installed')
                profiles[key+'_'+coordinate] = cells
                draw_profile(axis, cells, 'mean', key)
                axis.set(xlabel=coordinate, ylabel='mean variant / installed')
        axes[0].legend(fontsize=7)
        plots += save_plot(fig, output/(name+'_isolated_residual_and_l2'))
    if not data:
        matched = [r for r in rows if r['genpt'] is not None and r['genpt'] > 0]
        fig, axes = plt.subplots(2, 2, figsize=(11, 7))
        fig.suptitle(name+': native matched response; sparse bins explicit in JSON', fontsize=10)
        for rownum, (coordinate, edges) in enumerate([('genpt', PT_EDGES), ('geneta', np.linspace(-1.6, 1.6, 9))]):
            for key in ['rawpt', 'stored_jtpt', *variants]:
                cells = profile(matched, coordinate, edges, key, 'genpt')
                profiles[coordinate+'_'+key] = cells
                for col, statistic in enumerate(('mean', 'width')):
                    draw_profile(axes[rownum, col], cells, statistic, key)
                    axes[rownum, col].set(xlabel=coordinate, ylabel='response '+statistic)
        axes[0, 0].legend(fontsize=7)
        plots += save_plot(fig, output/(name+'_native_response'))
    return profiles, plots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('evidence', 'upstream', 'output'):
        parser.add_argument('--'+key, type=Path, required=True)
    args = parser.parse_args()
    manifest_path = args.evidence/'artifact_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    loaded = {name: load_sample(args.evidence, manifest, name, data) for name, (_, data) in SAMPLES.items()}
    chains = {name: payload_chains(system, data, args.upstream) for name, (system, data) in SAMPLES.items()}
    for variants in chains.values():
        for files in variants.values():
            for path in files:
                if not path.is_file():
                    raise ValueError('Missing retained payload: '+str(path))
    args.output.mkdir(parents=True, exist_ok=False)
    command = setup_library(args.output)
    summary = {'status': 'BOUNDED_RETAINED_LEDGER_KINEMATIC_DIAGNOSTIC_ONLY', 'limits': LIMITS,
               'manifest': {'path': str(manifest_path.resolve()), 'sha256': digest(manifest_path)},
               'compile_command': command, 'samples': [],
               'source_sha256': {str(p.relative_to(REPO)): digest(p) for p in [Path(__file__).resolve(), REPO/'validation/run3Jets/compare_calibrations.py', REPO/'src/JetCorrector.cxx', REPO/'src/JetCorrector.h']}}
    for name, (system, data) in SAMPLES.items():
        events, meta, provenance = loaded[name]
        correctors = {key: Corrector(files) for key, files in chains[name].items()}
        rows, counts = compare_events(events, system, data, correctors)
        ledger = args.output/(name+'_paired_jets.jsonl')
        ledger.write_text(''.join(json.dumps(row, allow_nan=False)+'\n' for row in rows))
        profiles, plots = plot_sample(rows, name, data, chains[name], args.output)
        summary['samples'].append({'label': name, 'system': system, 'data': data, 'counts': counts,
            'native_provenance': {k: meta[k] for k in ('input', 'file_uuid', 'expected_tree', 'event_tree', 'events_inspected', 'production_configuration_match', 'forest_metadata')},
            'retained_inputs': provenance, 'entry_range': [0, len(events)],
            'payloads': {key: [{'path': str(p.resolve()), 'sha256': digest(p)} for p in files] for key, files in chains[name].items()},
            'ledger': {'path': str(ledger.resolve()), 'sha256': digest(ledger)}, 'profiles': profiles, 'plots': plots})
        print(name, json.dumps(counts, sort_keys=True), flush=True)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
