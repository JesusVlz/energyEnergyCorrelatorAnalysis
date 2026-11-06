"""Build inputs must be explicit and bind the packaged native binary."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

import prepare_campaign as prepare


class BuildStateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.snapshot = self.folder/'snapshot'
        self.snapshot.mkdir()
        (self.snapshot/'eecAnalysis').write_bytes(b'native fixture')
        self.build = dict(status='BUILT_BATCH_RUNTIME', snapshot='snapshot',
                          binary_sha256=prepare.digest(self.snapshot/'eecAnalysis'),
                          cmssw_base='/cvmfs/fixture/CMSSW', scram_arch='el9_fixture')
        self.path = self.folder/'build.json'
        self.path.write_text(json.dumps(self.build))

    def test_relative_snapshot_resolves_beside_build_state(self):
        actual = prepare.load_build_state(self.path)
        self.assertEqual(actual['snapshot'], str(self.snapshot))
        self.assertEqual(actual['binary_sha256'], self.build['binary_sha256'])

    def test_changed_native_binary_is_rejected(self):
        (self.snapshot/'eecAnalysis').write_bytes(b'changed native fixture')
        with self.assertRaisesRegex(ValueError, 'missing or changed'):
            prepare.load_build_state(self.path)

    def test_absent_prerequisites_have_actionable_messages(self):
        with self.assertRaisesRegex(ValueError, '--build-state'):
            prepare.load_build_state(self.folder/'missing.json')
        with self.assertRaisesRegex(ValueError, '--inventory'):
            prepare.freeze(self.folder/'missing-inventory.json', 'fixture')

    def test_historical_dag_rejects_mismatched_runtime_before_outputs(self):
        args = SimpleNamespace(state_dir=self.folder/'new-state', build_state=self.path)
        with self.assertRaisesRegex(ValueError, 'DAG worker requires CMSSW_15_0_11'):
            prepare.prepare(args)
        self.assertFalse(args.state_dir.exists())

    def test_cli_accepts_explicit_state_and_fails_before_preparation(self):
        result = subprocess.run([sys.executable, str(Path(prepare.__file__)),
                                 '--state-dir', str(self.folder/'new-state'),
                                 '--build-state', str(self.folder/'missing.json')],
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('Missing build-state file', result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        self.assertFalse((self.folder/'new-state').exists())


if __name__ == '__main__':
    unittest.main()
