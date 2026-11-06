#!/usr/bin/env python3
"""Independently compare native pilot donor ledgers with bounded input-audit rows.

Checks exact 64-bit keys, source order, NH, untouched raw HF, fingerprint fold,
and global signal holdout flags. Final reco-veto decisions remain native.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def load(path):
    return [json.loads(line) for line in path.open()]


def verify(audit_dir, native_pool, signal_cap, pool_cap):
    manifest=json.loads((audit_dir/'manifest.json').read_text())
    signal=load(audit_dir/'signal.events.jsonl')[:signal_cap]
    if len(signal)!=signal_cap:
        raise ValueError('Audit signal range does not cover the requested native cap')
    ids={tuple(r['id']) for r in signal}
    backgrounds={r['background_key'] for r in signal if r['background_key']}
    offsets={}; offset=0
    for row in manifest['pool_files']:
        offsets[row['path']]=offset;offset+=row['source_entries']
    donors=load(audit_dir/'poolA.events.jsonl')+load(audit_dir/'poolB.events.jsonl')
    expected={tuple(r['id']):r for r in donors if offsets[r['source_file']]+r['source_entry']<pool_cap}
    if len(expected)!=pool_cap:
        raise ValueError('Audit pool range does not cover the requested native cap or has duplicate physical IDs')
    seen=set();errors=[];counts=Counter()
    def check(identity,field,actual,wanted):
        if actual!=wanted:errors.append({'id':list(identity),'field':field,'actual':actual,'expected':wanted})
    with native_pool.open() as stream:
        for row in csv.DictReader(stream):
            identity=tuple(int(row[k]) for k in ('run','lumi','event'))
            counts['rows']+=1
            if identity in seen: errors.append({'id':list(identity),'field':'duplicate_native_identity'})
            seen.add(identity)
            if identity not in expected:
                errors.append({'id':list(identity),'field':'identity_outside_audit_bound'})
                continue
            r=expected[identity]
            for field,wanted in [('donor_key',int(r['event_key'],16)),('nh',r['nh']),('hiBin',r['hiBin']),
                                 ('entry',offsets[r['source_file']]+r['source_entry']),
                                 ('background_key',int(r['background_key'],16) if r['background_key'] else 0)]:
                check(identity,field,int(row[field]),wanted)
            for field,name in [('vz','vz'),('hf','hf_raw'),('generator_weight','generator_weight')]:
                check(identity,field,float(row[field]),r[name])
            for field,wanted in [('pool',r['pool']),('background_signature',r['background_kinematics_signature'])]:
                check(identity,field,row[field],wanted)
            physical_excluded=identity in ids
            background_excluded=bool(r['background_key'] and r['background_key'] in backgrounds)
            for field,wanted in [('excluded_signal_physical_id',physical_excluded),('excluded_signal_background',background_excluded)]:
                counts[field]+=int(row[field]);check(identity,field,bool(int(row[field])),wanted)
            eligible=bool(int(row['eligible']));counts['eligible']+=eligible
            counts['pre_veto_source_selection_pass']+=r['pre_veto_selection_pass']
            if eligible:
                counts['eligible_pool_'+row['pool']]+=1
                if physical_excluded or background_excluded or not r['pre_veto_selection_pass']:
                    errors.append({'id':list(identity),'field':'invalid_eligible_event'})
            elif r['pre_veto_selection_pass'] and not physical_excluded and not background_excluded:
                counts['remaining_native_veto_or_selection_rejections']+=1
    check((), 'native_row_count',counts['rows'],pool_cap)
    check((), 'native_unique_identity_count',len(seen),pool_cap)
    return {'status':'PASS' if not errors else 'FAIL','scope':'Independent full bounded native donor ledger versus exact audit source identities and generator matching variables',
            'signal_holdout_cap':signal_cap,'pool_cap':pool_cap,'counts':dict(counts),'error_count':len(errors),'errors':errors[:100],
            'exact_fields':['physical FNV key','full run/lumi/event','global source entry','NH','hiBin','vz','raw hiHF_pf',
                            'raw generator weight','observed background key/signature','A/B partition','global signal physical/background holdout flags'],
            'limitations':'Native reconstructed jet veto not independently recomputed here; eligible entries must satisfy independently audited pre-veto source selection and all global identity exclusions.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audit_dir',type=Path)
    parser.add_argument('native_pool',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--signal-events',type=int,required=True)
    parser.add_argument('--pool-events',type=int,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Fresh verification output required')
    result=verify(args.audit_dir,args.native_pool,args.signal_events,args.pool_events)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['status']=='PASS' else 1)


if __name__=='__main__':main()
