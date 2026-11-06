"""Local worker regressions; no ROOT data, CRAB submission, or network access."""
import io
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest

WRAPPER = Path(__file__).with_name('compileAndRun.sh').resolve()


class WorkerTests(unittest.TestCase):
    def run_worker(self, mode='success', mixing='none'):
        temporary = tempfile.TemporaryDirectory(prefix='oo-worker-test-')
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        (directory / 'card.input').write_text('submitted card\n')
        (directory / 'FrameworkJobReport.xml').write_text('original framework report\n')
        with tarfile.open(directory / 'eecOO.tar.gz', 'w:gz') as archive:
            for name, data, permission in (
                ('card.input', b'obsolete archived card\n', 0o644),
                ('eecAnalysis', b'#!/bin/bash\ntouch "$CASE_DIRECTORY/stale_executable_ran"\n', 0o755),
            ):
                info = tarfile.TarInfo(name)
                info.size, info.mode = len(data), permission
                archive.addfile(info, io.BytesIO(data))
        with tarfile.open(directory / 'input_files.tar.gz', 'w:gz') as archive:
            data = b'["/store/test.root"]\n'
            info = tarfile.TarInfo('job_input_file_list_1.txt')
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
        bindir = directory / 'bin'
        bindir.mkdir()
        make = bindir / 'make'
        make.write_text('''#!/bin/bash
set -eu
if [[ $1 == -C ]]; then cd "$2"; shift 2; fi
if [[ $1 == clean ]]; then
  rm -f eecAnalysis
  exit 0
fi
if [[ $TEST_MODE == build_failure ]]; then exit 42; fi
cp "$CASE_DIRECTORY/fake_analyzer" eecAnalysis
chmod +x eecAnalysis
''')
        make.chmod(0o755)
        analyzer = directory / 'fake_analyzer'
        analyzer.write_text('''#!/bin/bash
set -eu
touch "$CASE_DIRECTORY/analyzer_ran"
cat "$2" > "$CASE_DIRECTORY/card_seen.txt"
printf '%s\\n' "$@" > "$CASE_DIRECTORY/args.txt"
if [[ $TEST_MODE == analyzer_failure ]]; then exit 17; fi
if [[ $TEST_MODE == missing_output ]]; then exit 0; fi
printf 'test output\\n' > "$3"
''')
        analyzer.chmod(0o755)
        environment = dict(os.environ, PATH=str(bindir)+os.pathsep+os.environ['PATH'],
                           CASE_DIRECTORY=str(directory), TEST_MODE=mode)
        result = subprocess.run(['bash', str(WRAPPER), '1', 'card=card.input',
                                 'output=result.root', 'location=2', 'mixinglist='+mixing],
                                cwd=directory, env=environment, capture_output=True, text=True)
        return directory, result

    def test_submitted_card_wins_and_output_is_staged(self):
        directory, result = self.run_worker()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((directory/'card_seen.txt').read_text(), 'submitted card\n')
        self.assertEqual((directory/'card.input').read_text(), 'submitted card\n')
        self.assertTrue((directory/'result.root').is_file())
        self.assertFalse((directory/'stale_executable_ran').exists())
        self.assertEqual((directory/'FrameworkJobReport.xml').read_text(), 'original framework report\n')
        self.assertEqual((directory/'args.txt').read_text().splitlines()[-1], '0')

    def test_numeric_minus_one_is_preserved(self):
        directory, result = self.run_worker(mixing='-1')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((directory/'args.txt').read_text().splitlines()[-1], '-1')

    def test_build_failure_never_runs_any_analyzer(self):
        directory, result = self.run_worker(mode='build_failure')
        self.assertEqual(result.returncode, 42)
        self.assertFalse((directory/'analyzer_ran').exists())
        self.assertFalse((directory/'stale_executable_ran').exists())
        self.assertFalse((directory/'result.root').exists())

    def test_analyzer_failure_is_propagated(self):
        directory, result = self.run_worker(mode='analyzer_failure')
        self.assertEqual(result.returncode, 17)
        self.assertFalse((directory/'result.root').exists())

    def test_missing_output_is_a_failure(self):
        directory, result = self.run_worker(mode='missing_output')
        self.assertEqual(result.returncode, 65)
        self.assertFalse((directory/'result.root').exists())


if __name__ == '__main__':
    unittest.main()
