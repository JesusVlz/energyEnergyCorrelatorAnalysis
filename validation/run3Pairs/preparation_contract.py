"""Small preparation contract shared by native exports and pair-table writers.

This records technical compatibility, never a physics-closure certificate.
"""
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
VERSION = 'run3-tracking-2026-11-v1'
PAIR_VERSION = 'run3-an751-v2'
VERTEX = {'Nominal': 3., 'Loose': 5., 'Tight': 2.}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def card_values(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split('#', 1)[0].split()
        if len(fields) >= 2:
            if fields[0] in result:
                raise ValueError('Repeated card key: ' + fields[0])
            result[fields[0]] = ' '.join(fields[1:])
    return result


def selection_from_card(path):
    values = card_values(path)
    error = float(values['MaxTrackPtRelativeError'])
    vertex = float(values['VertexMaxDistance'])
    inferred = next((name for name, limit in VERTEX.items() if math.isclose(vertex, limit, abs_tol=1e-6, rel_tol=0)), None)
    if not math.isclose(error, .1, abs_tol=1e-6, rel_tol=0) or inferred is None:
        raise ValueError('Invalid Run-3 cuts: require relative-error limit 0.1 and vertex limit 3/5/2')
    selection = values.get('Run3TrackSelection', inferred)
    if selection != inferred:
        raise ValueError('Run3TrackSelection disagrees with numeric cuts')
    return selection


def payload_path(system, selection):
    if selection not in VERTEX or system not in ('OO', 'pp'):
        raise ValueError('Unknown system or Run-3 track selection')
    if system == 'OO':
        return REPO / ('trackCorrectionTables/OO2025/Eff_OO_2025_Hijing_MB_Centrality_fromHihfpf_NoPU_3D_' + selection + '_Official_18Nov2025.root')
    return REPO / ('trackCorrectionTables/pp2024/Eff_ppref_2024_Pythia_QCD_pThat15_NopU_2D_' + selection + '_Official.root')


def export_contract(system, selection, payload):
    return {'version': VERSION, 'system': system, 'selection': selection,
            'payload_sha256': digest(payload), 'physical_charge': True,
            'coordinates': 'detector eta; source-event hiBin',
            'momentum_error_cut': 'pt > 10 => sigmaPt/pt <= 0.1; nonnegative finite error below',
            'vertex_limit_inclusive': VERTEX[selection],
            'oo_lookup_pt_cap_GeV': 6.1 if system == 'OO' else None}


def validate_contract(contract):
    if not isinstance(contract, dict) or contract.get('version') != VERSION:
        raise ValueError('Stale Run-3 tracking input: regenerate with the corrected native selection')
    if contract.get('selection') not in VERTEX or contract.get('system') not in ('OO', 'pp'):
        raise ValueError('Missing system/selection in tracking contract')
    checksum = contract.get('payload_sha256', '')
    if len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum):
        raise ValueError('Tracking contract requires the actual single-track payload SHA256')
    if contract.get('physical_charge') is not True:
        raise ValueError('Tracking input does not declare physical reconstructed charge')
    return contract


def check_export_metadata(paths):
    """Reject stale native exports and mixtures before computing any pair map."""
    by_system = {}
    for path in paths:
        path = Path(path)
        meta = json.loads(Path(str(path) + '.metadata.json').read_text())
        contract = validate_contract(meta.get('tracking_contract'))
        if meta.get('output_sha256') != digest(path):
            raise ValueError('Export checksum differs from receipt: ' + str(path))
        system = contract['system']
        if system in by_system and by_system[system] != contract:
            raise ValueError('Mixed track selections/payloads in pair derivation for ' + system)
        by_system[system] = contract
    if len({value['selection'] for value in by_system.values()}) != 1:
        raise ValueError('Keep Nominal/Loose/Tight derivations separate')
    return by_system


def write_contract(ROOT, result, rootfile, system):
    contract = result.get('tracking_contracts', {}).get(system)
    if contract is None:
        # Historical/synthetic test results may still be drawn but cannot be
        # accepted by the repaired analyzer as a new correction payload.
        ROOT.TNamed('Run3TrackSelectionVersion', 'UNVERIFIED_HISTORICAL_INPUT').Write()
        return
    validate_contract(contract)
    ROOT.TNamed('Run3TrackSelectionVersion', VERSION).Write()
    ROOT.TNamed('Run3TrackSelection', contract['selection']).Write()
    ROOT.TNamed('Run3PairCollisionSystem', system).Write()
    ROOT.TNamed('Run3SingleTrackPayloadSHA256', contract['payload_sha256']).Write()
