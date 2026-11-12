#!/usr/bin/env python3
"""Plan or submit the R4 MC inputs needed for the OO/ppRef EEC workflow.

The OO measurement cards include the two same-event cone components. The
ppRef MC measurement has no background cones. The existing HIJING pool A
requires event-level reservations and belongs to the separate strict GenGen
closure workflow; ordinary mixed-cone MC is not in this CRAB campaign.
Pair preparation produces inputs for a later derived correction;
it does not enable an unqualified pair correction in the measurement cards.
"""

import argparse
import subprocess
import sys
from pathlib import Path

from run3_crab import CRAB_DIR, card_for, files_for


REPO = CRAB_DIR.parent


def cards(sample):
    system, scopes, pair_prefix = (
        ('OO', ('integrated', 'split'), 'OO') if sample == 'oo-mc'
        else ('ppRef', ('inclusive',), 'pp')
    )
    core = []
    for scope in scopes:
        stem = f'{system}_MC_{scope}'
        cone_mode = 'sameEventCones' if sample == 'oo-mc' else 'none'
        core += [f'cards/R4/measurement/{stem}_{cone_mode}.input']
        core += [f'cards/R4/mc-closure/{stem}_{mode}.input'
                 for mode in ('RecoGen', 'GenReco', 'GenGen')]
        core += [f'cards/R4/response/{stem}_truthGuard5.input',
                 f'cards/R4/jet-qa/{stem}_jetResponse.input',
                 f'cards/R4/particle-response/{stem}_trackPairResolution.input']
    pairs = [f'cards/R4/pairs/{pair_prefix}_{selection}_{mode}.input'
             for selection in ('Nominal', 'Loose', 'Tight')
             for mode in ('RecoReco', 'RecoGen')]
    return core, pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', choices=('all', 'oo-mc', 'ppref-mc'), default='all')
    parser.add_argument('--stage', choices=('all', 'core', 'pairs'), default='all')
    parser.add_argument('--tag', default='R4v1')
    parser.add_argument('--pilot', type=int, metavar='N')
    parser.add_argument('--submit', action='store_true')
    args = parser.parse_args()
    if args.pilot is not None and not 1 <= args.pilot <= 100:
        parser.error('--pilot must select 1..100 files')
    samples = ('oo-mc', 'ppref-mc') if args.sample == 'all' else (args.sample,)
    planned = []
    for sample in samples:
        core, pairs = cards(sample)
        selected = ((core if args.stage in ('all', 'core') else [])
                    + (pairs if args.stage in ('all', 'pairs') else []))
        for choice in selected:
            card, values = card_for(sample, choice)
            files_for(sample, 1, 1, values)
            planned.append((sample, card))
    if args.submit and not (CRAB_DIR / 'eecOO.tar.gz').is_file():
        parser.error('Build crab/eecOO.tar.gz before submitting')
    print(f'{len(planned)} MC tasks to {"submit" if args.submit else "plan"}', flush=True)
    for sample, card in planned:
        command = [sys.executable, str(CRAB_DIR / 'submitRun3.py'), sample,
                   '--card', str(card.relative_to(REPO)), '--tag', args.tag]
        if args.pilot is not None:
            command += ['--pilot', str(args.pilot)]
        if args.submit:
            command.append('--submit')
        print(' '.join(command), flush=True)
        subprocess.run(command, cwd=REPO, check=True)


if __name__ == '__main__':
    main()
