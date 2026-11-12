#!/usr/bin/env python3
"""Plan Run-3 EEC tasks and, only with --submit, send them to CRAB."""

import argparse
import os
import subprocess

from run3_crab import CRAB_DIR, card_for, files_for, task_name

CONFIGS = {
    'oo-data': 'crabEECOO.py', 'ppref-data': 'crabEECppref5v36.py',
    'oo-mc': 'crabEECOOMC.py', 'ppref-mc': 'crabEECpprefMC5v36.py',
}


def tasks(sample):
    if sample == 'oo-data':
        return [(part, chunk) for part in range(1, 11) for chunk in (1, 2)]
    if sample == 'ppref-data':
        return [(part, 1) for part in range(1, 5)]
    return [(1, 1)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample', choices=CONFIGS)
    parser.add_argument('--card', help='Runnable repository-relative cards/R4/*.input card')
    parser.add_argument('--tag', default='R4v1', help='Unique campaign revision, e.g. R4v2')
    parser.add_argument('--parts', help='Comma-separated data part numbers; omit for all')
    parser.add_argument('--pilot', type=int, metavar='N', help='Use N files from one task as a separate pilot')
    parser.add_argument('--submit', action='store_true', help='Actually invoke crab submit')
    args = parser.parse_args()
    if args.pilot is not None and not 1 <= args.pilot <= 100:
        parser.error('--pilot must select 1..100 files')
    if args.pilot is not None and args.parts:
        parser.error('Select --pilot or --parts, not both')
    # The driver owns the override; do not inherit one from a previous shell.
    os.environ.pop('EEC_INPUT_LIST', None)
    card, values = card_for(args.sample, args.card)
    chosen = tasks(args.sample)[:1] if args.pilot is not None else tasks(args.sample)
    if args.parts:
        requested = {int(x) for x in args.parts.split(',')}
        if not requested or requested - {part for part, _ in chosen}:
            parser.error('Unknown part number in --parts')
        chosen = [(part, chunk) for part, chunk in chosen if part in requested]
    names = []
    all_files = set()
    tag = args.tag + (f'pilot{args.pilot}' if args.pilot is not None else '')
    for part, chunk in chosen:
        files = files_for(args.sample, part, chunk, values)
        if args.pilot is not None:
            files = files[:args.pilot]
        if all_files.intersection(files):
            raise ValueError(f'Duplicate signal file across tasks at part {part} chunk {chunk}')
        all_files.update(files)
        name = task_name(args.sample, card, part, chunk, tag)
        names.append((name, part, chunk, files))
    if args.submit and not (CRAB_DIR / 'eecOO.tar.gz').is_file():
        raise ValueError('Build crab/eecOO.tar.gz first: python3 crab/makeOOAnalysisTar.py')
    existing = [str(CRAB_DIR / 'crab_projects' / 'run3' / ('crab_' + name))
                for name, _, _, _ in names
                if (CRAB_DIR / 'crab_projects' / 'run3' / ('crab_' + name)).exists()]
    if existing:
        raise ValueError('Existing CRAB project(s); choose a new tag or review them:\n'
                         + '\n'.join(existing))
    for name, part, chunk, files in names:
        env = os.environ.copy()
        env.update({'EEC_CARD': str(card), 'EEC_RUN_TAG': tag})
        if args.pilot is not None:
            pilot_list = CRAB_DIR / 'inputLists' / (name + '.txt')
            if args.submit:
                pilot_list.parent.mkdir(exist_ok=True)
                pilot_list.write_text('\n'.join(files) + '\n')
            env['EEC_INPUT_LIST'] = str(pilot_list)
        if args.sample == 'oo-data':
            env.update({'OO_CRAB_PART': str(part), 'OO_CRAB_CHUNK': str(chunk)})
        elif args.sample == 'ppref-data':
            env.update({'PPREF_CRAB_PART': str(part), 'PPREF_CRAB_CHUNK': '1'})
        count = len(files)
        print(f'{name}: {count} files, {count if args.sample.endswith("data") else (count+1)//2} jobs', flush=True)
        if args.submit:
            subprocess.run(['crab', 'submit', '-c', CONFIGS[args.sample]],
                           cwd=CRAB_DIR, env=env, check=True)
    print(f'{len(names)} tasks, {len(all_files)} distinct files; '
          + ('submitted' if args.submit else 'plan only (add --submit to submit)'))


if __name__ == '__main__':
    main()
