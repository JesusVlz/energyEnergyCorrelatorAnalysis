#!/usr/bin/env python3
"""Independently validate bounded native OO closure ledgers against audited inputs."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import re

SIGNAL_NAMESPACE = 'OO2025_embedded_existing_owner_pilot'
POOL_NAMESPACE = 'OO2025_HIJING_existing_owner_pool_pilot'
BACKGROUND_NAMESPACE = 'OO2025_HIJING_stored_kinematics_v1'
MASK = (1 << 64)-1


def require(condition, message):
    if not condition:
        raise ValueError(message)


def uint64(value):
    require(isinstance(value, (str, int)) and re.fullmatch(r'[0-9]+', str(value)) is not None,
            'Identity/key must be a decimal integer without float conversion: '+str(value))
    number = int(value)
    require(0 <= number <= MASK, 'Identity exceeds uint64')
    return number


def event_key(namespace, identity):
    key = 14695981039346656037
    payload = namespace.encode()+b'\0'+b''.join(uint64(v).to_bytes(8, 'little') for v in identity)
    for byte in payload:
        key = ((key ^ byte)*1099511628211) & MASK
    return key


def fingerprint_key(signature):
    count, total, xor = signature.split(':')
    return event_key(BACKGROUND_NAMESPACE, [int(count), int(total, 16), int(xor, 16)]) if int(count) else 0


def identity(row, prefix=''):
    return tuple(uint64(row[prefix+name]) for name in ('run', 'lumi', 'event'))


def finite(value):
    result = float(value)
    require(math.isfinite(result), 'Nonfinite ledger value: '+str(value))
    return result


def rows(path):
    with Path(path).open(newline='') as source:
        reader = csv.DictReader(source)
        require(reader.fieldnames is not None, 'Missing CSV header: '+str(path))
        for row in reader:
            require(None not in row and None not in row.values(), 'CSV field count mismatch: '+str(path))
            yield row


def jsonlines(path):
    with Path(path).open() as source:
        for line in source:
            if line.strip():
                yield json.loads(line)


def same(number, expected, description):
    require(math.isclose(finite(number), finite(expected), rel_tol=1e-13, abs_tol=1e-13), description)


def validate(prefix, audit_dir, signal_events=10000, pool_events=100000, seed=20261030):
    prefix, audit_dir = str(prefix), Path(audit_dir)
    manifest = json.loads((audit_dir/'manifest.json').read_text())
    require(signal_events > 0 and pool_events > 0, 'Explicit positive global caps required')
    signal_rows = list(jsonlines(audit_dir/'signal.events.jsonl'))[:signal_events]
    require(len(signal_rows) == min(signal_events, manifest['signal']['source_entries']),
            'Authoritative signal audit is shorter than the requested cap')
    signals = {tuple(row['id']): row for row in signal_rows}
    require(len(signals) == len(signal_rows), 'Duplicate audited signal identities')
    signal_bg = {fingerprint_key(row['background_kinematics_signature']) for row in signal_rows
                 if int(row['background_kinematics_signature'].split(':')[0])}
    offsets, offset = {}, 0
    for file in manifest['pool_files']:
        offsets[file['path']] = offset
        offset += file['source_entries']
    audited_pool = {}
    for group in ('A', 'B'):
        for row in jsonlines(audit_dir/('pool'+group+'.events.jsonl')):
            entry = offsets[row['source_file']]+row['source_entry']
            if entry < pool_events:
                audited_pool[tuple(row['id'])] = (entry, row)
    expected_pool = min(pool_events, offset)
    require(len(audited_pool) == expected_pool, 'Authoritative pool audit is shorter than requested global cap')

    donors, exclusions, occupancy = {}, Counter(), Counter()
    for row in rows(prefix+'.pool.csv'):
        physical = identity(row)
        require(physical in audited_pool, 'Native donor outside bounded audited pool')
        entry, audit = audited_pool[physical]
        key = uint64(row['donor_key'])
        require(key == event_key(POOL_NAMESPACE, physical), 'Donor uint64 key mismatch')
        require(key not in donors, 'Duplicate native donor physical identity/key')
        require(int(row['entry']) == entry, 'Donor source/global entry mismatch')
        same(row['vz'], audit['vz'], 'Donor vz differs from source')
        same(row['hf'], audit['hf_raw'], 'HF OFF changed the raw stored coordinate')
        same(row['generator_weight'], audit['generator_weight'], 'Donor weight differs from source')
        require(finite(row['generator_weight']) == 1., 'Nonunity donor sampling requires an explicit weighting study')
        require(int(row['nh']) == audit['nh'], 'Donor matching NH differs from fixed source definition')
        require(int(row['hiBin']) == audit['hiBin'], 'Donor physical hiBin changed')
        require(row['background_signature'] == audit['background_kinematics_signature'], 'Donor source fingerprint mismatch')
        bg = fingerprint_key(row['background_signature'])
        require(uint64(row['background_key']) == bg, 'Donor background key mismatch')
        group = 'B' if ((bg or key) ^ seed) % 2 else 'A'
        require(row['pool'] == group == audit['pool'], 'A/B partition differs from audited background identity')
        raw_excluded, bg_excluded = physical in signals, bg != 0 and bg in signal_bg
        require(int(row['excluded_signal_physical_id']) == raw_excluded, 'Global signal raw-ID exclusion mismatch')
        require(int(row['excluded_signal_background']) == bg_excluded, 'Global signal background exclusion mismatch')
        eligible = int(row['eligible']) == 1
        require(not eligible or not (raw_excluded or bg_excluded), 'Signal/donor leakage: excluded donor remains eligible')
        require(not eligible or audit['pre_veto_selection_pass'], 'Eligible donor failed audited source event selection')
        exclusions['physical_id'] += raw_excluded
        exclusions['observed_background'] += bg_excluded
        exclusions['union'] += raw_excluded or bg_excluded
        occupancy[group+'_scanned'] += 1
        occupancy[group+'_eligible'] += eligible
        donors[key] = row
    require(len(donors) == expected_pool, 'Native donor scan does not equal bounded source cap')
    by_nh = defaultdict(list)
    for key, row in donors.items():
        if int(row['eligible']):
            by_nh[int(row['nh'])].append((key, row))

    assignments, uses, failures = {}, Counter(), Counter()
    for row in rows(prefix+'.partners.csv'):
        physical = identity(row)
        require(physical in signals, 'Native signal outside capped audit')
        audit = signals[physical]
        key = uint64(row['signal_key'])
        require(key == event_key(SIGNAL_NAMESPACE, physical), 'Signal uint64 key mismatch')
        require(key not in assignments, 'More than one assignment per physical signal event')
        require(int(row['entry']) == audit['source_entry'], 'Signal entry differs from audited source')
        require(int(row['nh']) == audit['nh'] and int(row['hiBin']) == audit['hiBin'], 'Signal NH/hiBin mismatch')
        same(row['vz'], audit['vz'], 'Signal vz differs from source')
        same(row['raw_generator_weight'], audit['generator_weight'], 'Raw signal weight changed')
        finite(row['event_weight'])
        require(row['signal_background_signature'] == audit['background_kinematics_signature'], 'Signal fingerprint differs from source')
        bg = fingerprint_key(row['signal_background_signature'])
        require(uint64(row['signal_background_key']) == bg, 'Signal background key mismatch')
        require(int(row['inspected']) == expected_pool, 'Partner search was not exactly one finite capped scan')
        require(row['status'] in ('matched', 'failed'), 'Unknown assignment status')
        available = Counter()
        tolerance = max(1., .01*int(row['nh']))
        for nh in range(max(0, math.ceil(int(row['nh'])-tolerance)), math.floor(int(row['nh'])+tolerance)+1):
            for _, donor in by_nh[nh]:
                if abs(finite(donor['vz'])-finite(row['vz'])) <= .5:
                    available[donor['pool']] += 1
        require((row['status'] == 'matched') == bool(available['A'] and available['B']),
                'Native matching pass/failure disagrees with independent fixed-tolerance feasibility')
        if row['status'] == 'matched':
            keys = []
            for partner in (1, 2):
                base = 'donor'+str(partner)+'_'
                donor_key = uint64(row[base+'key'])
                require(donor_key in donors, 'Partner outside bounded native pool')
                donor = donors[donor_key]
                require(int(donor['eligible']) == 1, 'Partner is not eligible')
                require(identity(row, base) == identity(donor), 'Partner physical ID mismatch')
                require(int(row[base+'entry']) == int(donor['entry']), 'Partner source/global entry mismatch')
                same(row[base+'vz'], donor['vz'], 'Partner vz mismatch')
                require(int(row[base+'nh']) == int(donor['nh']), 'Partner NH mismatch')
                same(row[base+'weight'], donor['generator_weight'], 'Partner source weight mismatch')
                require(row[base+'background_key'] == donor['background_key'], 'Partner background key mismatch')
                require(donor['pool'] == ('A' if partner == 1 else 'B'), 'Partners do not use disjoint A/B pools')
                require(abs(finite(row[base+'vz'])-finite(row['vz'])) <= .5, 'Partner violates fixed dz')
                require(abs(int(row[base+'nh'])-int(row['nh'])) <= tolerance, 'Partner violates fixed NH tolerance')
                keys.append(donor_key)
                uses[donor_key] += 1
            require(keys[0] != keys[1], 'M1 and M2 identical')
            require(not (row['donor1_background_key'] != '0' and row['donor1_background_key'] == row['donor2_background_key']), 'M1/M2 share observed background')
        else:
            for partner in (1, 2):
                base = 'donor'+str(partner)+'_'
                require(all(row[base+field] == '' for field in ('key','run','lumi','event','entry','vz','nh','weight')),
                        'Failed event retained a partial estimator')
            failures[row['reason']] += 1
        assignments[key] = row
    for key, row in donors.items():
        require(int(row['uses']) == uses[key], 'Donor reuse counter differs from event assignments')

    jets, exposure, jet_counts = {}, Counter(), Counter()
    for row in rows(prefix+'.jets.csv'):
        key = uint64(row['signal_key'])
        require(key in assignments, 'Jet has no event-level partner assignment')
        require(row['status'] == assignments[key]['status'], 'Jet did not reuse its event assignment')
        index = (key, int(row['jet_index']))
        require(index not in jets, 'Duplicate jet exposure')
        for field in ('jet_pt','jet_eta','jet_phi','centrality','event_weight'):
            finite(row[field])
        same(row['event_weight'], assignments[key]['event_weight'], 'Jet/exposure event weights differ')
        require(abs(finite(row['jet_phi'])) <= math.pi+1e-6, 'Jet phi outside circular source range')
        jets[index] = row
        jet_counts[row['status']] += 1
        exposure[row['status']] += finite(row['event_weight'])

    contributions, powers, categories, total_value = 0, Counter(), Counter(), Counter()
    for row in rows(prefix+'.contributions.csv'):
        index = (uint64(row['signal_key']), int(row['jet_index']))
        require(index in jets and jets[index]['status'] == 'matched', 'Contribution outside common matched jet exposure')
        same(row['jet_pt'], jets[index]['jet_pt'], 'Contribution changed fixed jet pT')
        same(row['centrality'], jets[index]['centrality'], 'Contribution changed signal centrality')
        power = int(float(row['power']))
        require(power in (1,2) and float(row['power']) == power, 'Unexpected energy power')
        category, origin = int(row['pairing']), int(row['subevent'])
        require(category in (0,3,4,5,6,7,8,9,10,11,12), 'Unexpected pilot pairing category')
        require(0 <= origin <= 3, 'Invalid origin or aggregate PH counted a second time')
        if category in (4,6,7):
            require(origin == 3, 'Pure HIJING donor pair has non-background origin')
        if category in (3,5):
            require(origin in (1,3), 'Pure HIJING donor was incorrectly identified as PYTHIA')
        require(finite(row['delta_r_high']) > finite(row['delta_r_low']), 'Invalid delta-R bin bounds')
        require(int(row['delta_r_bin']) >= 0 and int(row['track_pt_bin']) >= 0, 'Negative ROOT bin index')
        value = finite(row['value'])
        contributions += 1
        powers[str(power)] += 1
        categories[str(category)] += 1
        total_value[str(power)] += value
    require(not contributions or set(powers) == {'1','2'}, 'One requested power missing from nonempty paired output')
    return {'status': 'PASS', 'prefix': prefix, 'audit_directory': str(audit_dir), 'seed': seed,
            'global_caps': {'signal': signal_events, 'pool': pool_events},
            'audited_signal_events': len(signals), 'native_pool_events': len(donors),
            'selected_signal_events': len(assignments),
            'matched_signal_events': sum(row['status']=='matched' for row in assignments.values()),
            'matching_failures': dict(failures), 'pool_occupancy': dict(occupancy),
            'global_holdout_exclusions': dict(exclusions), 'jet_counts': dict(jet_counts),
            'weighted_jet_exposure': dict(exposure), 'unique_assigned_partners': len(uses),
            'partner_assignments': sum(uses.values()), 'maximum_partner_reuse': max(uses.values(), default=0),
            'contribution_rows': contributions, 'rows_by_power': dict(powers), 'rows_by_pairing': dict(categories),
            'contribution_sum_by_power': dict(total_value),
            'checks': ['fixed global caps', 'full uint64 physical IDs/keys', 'global signal/donor identity and observed-background exclusions',
                       'raw HF exact identity', 'fixed NH independent of EEC cuts', 'background-based A/B partition',
                       'finite eligible-only exact-tolerance matching', 'no partial failures', 'shared event assignment across jets/powers',
                       'common matched jet exposure', 'pure HIJING classification', 'event-level reuse reconstruction', 'finite weights and unity donors'],
            'independence': 'Unknown production lineage remains a limitation. Observed identity/fingerprint disjointness is verified within bounded scanned ranges; partner reuse is not independent exposure.',
            'uncertainty': 'These checks do not establish final covariance or physics closure.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--audit-dir', type=Path, required=True)
    parser.add_argument('--signal-events', type=int, default=10000)
    parser.add_argument('--pool-events', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=20261030)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.prefix, args.audit_dir, args.signal_events, args.pool_events, args.seed)
    with args.output.open('x') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
