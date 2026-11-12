#!/usr/bin/env python3
"""List the actual OO and ppRef MC forests mounted on lxplus EOS."""

from pathlib import Path

from run3_crab import CRAB_DIR, EOS_PREFIXES


def main():
    output_dir = CRAB_DIR / 'inputLists'
    output_dir.mkdir(exist_ok=True)
    for sample, pattern, name in (
        ('oo-mc', 'unembedded_L1HLT_objects_MC_*.root', 'OOMC_EOS.txt'),
        ('ppref-mc', 'HiForestMiniAOD_MC_*.root', 'ppRefMC_EOS.txt'),
    ):
        prefix = EOS_PREFIXES[sample]
        # The EOS user XRootD namespace has a double slash before /eos/.
        directory = Path('/eos/' + prefix.split('//eos/', 1)[1])
        files = sorted(directory.glob(pattern), key=lambda p: (int(p.stem.rsplit('_', 1)[-1]), p.name))
        if not files:
            raise SystemExit(f'No {pattern} files found in {directory}; check the EOS mount')
        paths = [prefix + p.name for p in files]
        if len(set(paths)) != len(paths):
            raise SystemExit(f'Duplicate MC input files in {directory}')
        destination = output_dir / name
        destination.write_text('\n'.join(paths) + '\n')
        print(f'{sample}: {len(paths)} files -> {destination}')
        print('  first:', paths[0])
        print('  last: ', paths[-1])


if __name__ == '__main__':
    main()
