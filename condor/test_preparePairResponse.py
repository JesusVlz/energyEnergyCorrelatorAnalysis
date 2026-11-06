"""Preparation consumes explicit runtime inputs and optional campaign references."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import preparePairResponse as prepare


class PreparationTests(unittest.TestCase):
    def test_submit_description_uses_supplied_runtime(self):
        build = {'cmssw_base': '/cvmfs/fixture/CMSSW', 'scram_arch': 'el9_fixture'}
        production = prepare.submit_description(Path('/tmp/fixture'), 'a'*64, 'OOMC', 'OO', 'tag', build)
        merge = prepare.submit_description(Path('/tmp/fixture'), 'a'*64, 'ppMC', 'pp', 'tag', build, True, 'pairs')
        for text in (production, merge):
            self.assertIn('/cvmfs/fixture/CMSSW el9_fixture '+'a'*64, text)
        self.assertIn('queue mode,input_url', production)
        self.assertIn('queue 1', merge)

    def test_cli_documents_prerequisites_and_optional_references(self):
        result = subprocess.run([sys.executable, str(Path(prepare.__file__)), '--help'],
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for option in ('--build-state', '--state-dir', '--inventory', '--reference-oomc', '--reference-ppmc'):
            self.assertIn(option, result.stdout)

    def test_missing_build_state_is_actionable_before_any_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            result = subprocess.run([sys.executable, str(Path(prepare.__file__)),
                                     '--state-dir', str(folder/'state'),
                                     '--build-state', str(folder/'missing.json'),
                                     '--inventory', str(folder/'inventory.json'),
                                     '--condor-root', str(folder/'condor')],
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('Missing build-state file', result.stderr)
            self.assertNotIn('Traceback', result.stderr)
            self.assertFalse((folder/'state').exists())
            self.assertFalse((folder/'condor').exists())


if __name__ == '__main__':
    unittest.main()
