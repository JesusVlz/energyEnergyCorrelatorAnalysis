#!/usr/bin/env python3
"""Check/update active Run-3 cards. Historical evidence and Run-2 are immutable."""
import argparse
import json
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
OBSOLETE = {'CalorimeterSignalLimitPt', 'HighPtEtFraction', 'Chi2QualityCut', 'MinimumTrackHits'}
HISTORICAL = {'archive', 'evidence', 'results', 'pending'}

def fields(line):
    return line.split('#', 1)[0].split()

def inspect(write=False):
    checked, changed = [], []
    paths = sorted(set(REPO.rglob('*.input')) | set(REPO.rglob('*.input.in')))
    for path in paths:
        relative = path.relative_to(REPO)
        if any(part in HISTORICAL or part == '.git' for part in relative.parts):
            continue
        lines = path.read_text().splitlines(keepends=True)
        entries = [fields(line) for line in lines]
        kinds = [int(float(v[1])) for v in entries if len(v)>1 and v[0]=='DataType']
        if kinds not in ([9], [10], [11], [12]):
            continue
        checked.append(str(relative))
        cleaned = [line for line in lines if not fields(line) or fields(line)[0] not in OBSOLETE]
        if lines != cleaned:
            changed.append(str(relative))
            if write:
                path.write_text(''.join(cleaned))
    if changed and not write:
        raise ValueError('Obsolete Run-3 card keys remain: ' + ', '.join(changed))
    return {'run3_cards_checked':len(checked), 'updated':changed,
            'preserved':'Run-2 and archived/evidence/results/pending records'}

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write',action='store_true')
    print(json.dumps(inspect(p.parse_args().write),indent=2))
