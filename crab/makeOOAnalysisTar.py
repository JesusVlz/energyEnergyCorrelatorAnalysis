#!/usr/bin/env python3
"""Build the CRAB worker archive from current sources and calibration payloads.

The archive excludes configuration cards, executables, and object files. It
records the exact packaged contents in worker_manifest.json, including local
source changes that are not represented by the Git commit alone.
"""

import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile


def main():
    crab = Path(__file__).resolve().parent
    repository = crab.parent
    output = crab / 'eecOO.tar.gz'
    sources = [repository / 'Makefile', repository / 'eecAnalysis.cxx']
    sources += sorted(path for path in (repository / 'src').rglob('*')
                      if path.is_file() and path.suffix in ('.h', '.cxx', '.cc', '.hh', '.hpp'))
    for directory in ('jetEnergyCorrections', 'jetVetoMaps', 'trackCorrectionTables', 'mixingFileList'):
        sources += sorted(path for path in (repository / directory).rglob('*')
                          if path.is_file() and not any(part.startswith('.') for part in path.relative_to(repository).parts)
                          and path.suffix not in ('.o', '.so', '.pyc'))
    required = (
        'src/Run3JetSelection.cxx', 'src/Run3JetSelectionInputs.h',
        'src/Run3JetBranchValidation.h', 'src/Run3JetGeometry.h',
        'src/OOMixingPilot.h', 'src/OOEventIdentity.h',
        'jetEnergyCorrections/L2Relative_AK4PF_OO_V1.txt',
        'jetVetoMaps/OO2025/Summer24Prompt25_RunCDEFG.root',
        'mixingFileList/OO2025_MB_pool_pilot.txt',
    )
    names = {path.relative_to(repository).as_posix() for path in sources}
    missing = set(required) - names
    if missing:
        raise RuntimeError('Missing worker payloads: ' + ', '.join(sorted(missing)))
    revision = subprocess.check_output(['git', '-C', str(repository), 'rev-parse', 'HEAD'], text=True).strip()
    status = subprocess.check_output(
        ['git', '-C', str(repository), 'status', '--porcelain', '--', 'Makefile', 'eecAnalysis.cxx',
         'src', 'jetEnergyCorrections', 'jetVetoMaps', 'trackCorrectionTables', 'mixingFileList'], text=True)
    manifest = {'git_commit': revision, 'working_tree_dirty': bool(status.strip()), 'files_sha256': {}}
    with tempfile.NamedTemporaryFile(dir=crab, prefix='.eecOO.', suffix='.tar.gz', delete=False) as tmp:
        temporary = Path(tmp.name)
    try:
        with tarfile.open(temporary, 'w:gz') as archive:
            for path in sources:
                name = path.relative_to(repository).as_posix()
                data = path.read_bytes()
                if name == 'eecAnalysis.cxx':
                    data = data.replace(b'GITHASHHERE', revision.encode())
                manifest['files_sha256'][name] = hashlib.sha256(data).hexdigest()
                info = tarfile.TarInfo(name)
                info.size = len(data)
                info.mode = 0o644
                info.mtime = int(path.stat().st_mtime)
                archive.addfile(info, io.BytesIO(data))
            data = (json.dumps(manifest, indent=2, sort_keys=True) + '\n').encode()
            info = tarfile.TarInfo('worker_manifest.json')
            info.size = len(data)
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(data))
        if temporary.stat().st_size > 100 * 1024 * 1024:
            raise RuntimeError('Worker archive exceeds 100 MiB; review the CRAB sandbox size')
        os.chmod(temporary, 0o644)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f'Created {output.name}: {output.stat().st_size:,} bytes, {len(sources)} payload files')
    print('SHA256:', hashlib.sha256(output.read_bytes()).hexdigest())
    print('Git commit:', revision, '(local changes included)' if manifest['working_tree_dirty'] else '')
    print('No configuration cards, executable, or object files included.')


if __name__ == '__main__':
    main()
