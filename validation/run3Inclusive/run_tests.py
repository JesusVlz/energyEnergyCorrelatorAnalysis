#!/usr/bin/env python3
"""Build current sources in a temporary directory and run synthetic inclusive checks.

This deliberately avoids repository .o files, which may belong to older source revisions.
No production forests, submitted jobs, or physics-closure claims are involved.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]

def main():
    flags = shlex.split(subprocess.check_output(["root-config", "--cflags"], text=True))
    libraries = shlex.split(subprocess.check_output(["root-config", "--libs"], text=True))
    with tempfile.TemporaryDirectory(prefix="oo-inclusive-tests-") as tmp:
        build = Path(tmp)
        sources = sorted((ROOT / "src").glob("*.cxx"))
        def compile_source(source):
            target = build / (source.stem + ".o")
            subprocess.run(["g++", *flags, "-I", str(ROOT), "-c", str(source), "-o", str(target)], check=True)
            return str(target)
        with ThreadPoolExecutor(max_workers=4) as workers:
            objects = list(workers.map(compile_source, sources))
        executable = build / "checkInclusiveDiagnostics"
        subprocess.run(["g++", *flags, "-fno-access-control", "-I", str(ROOT),
                        str(ROOT / "validation/run3Inclusive/checkInclusiveDiagnostics.cxx"),
                        *objects, *libraries, "-lEG", "-o", str(executable)], check=True)
        subprocess.run([str(executable)], cwd=ROOT, check=True)

if __name__ == "__main__":
    main()
