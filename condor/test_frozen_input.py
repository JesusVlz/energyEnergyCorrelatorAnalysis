"""Validate frozen metadata across local and XRootD staging."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import frozen_input


class FrozenInputTest(unittest.TestCase):
    def test_xrootd_timestamp_is_utc_and_compares_at_protocol_precision(self):
        url = 'root://eos.example//eos/user/sample.root'
        expected = {'size_bytes': 42, 'mtime_ns': 1792735430485461000}
        result = SimpleNamespace(stdout='Path: /eos/user/sample.root\nSize: 42\nMTime: 2026-10-23 06:03:50\n')
        with mock.patch.object(frozen_input.subprocess, 'run', return_value=result) as command:
            actual = frozen_input.check_source_identity(url, expected)
        self.assertEqual(actual, {'size_bytes': 42, 'mtime_seconds': 1792735430})
        self.assertEqual(command.call_args.kwargs['env']['TZ'], 'UTC')
        self.assertEqual(command.call_args.args[0], ['xrdfs', 'root://eos.example', 'stat', '/eos/user/sample.root'])
        with mock.patch.object(frozen_input.subprocess, 'run', return_value=result):
            with self.assertRaisesRegex(ValueError, 'modification time changed'):
                frozen_input.check_source_identity(url, {**expected, 'mtime_ns': expected['mtime_ns'] - 1_000_000_000})

    def test_remote_change_during_staging_rejected_and_malformed_metadata_fails_closed(self):
        url = 'root://eos.example//eos/user/sample.root'
        expected = {'size_bytes': 3, 'mtime_ns': 10_000_000_000}
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / 'local.root'
            def stage():
                destination.write_bytes(b'abc')
            with mock.patch.object(frozen_input, 'source_metadata', side_effect=[
                    {'size_bytes': 3, 'mtime_seconds': 10}, {'size_bytes': 3, 'mtime_seconds': 11}]):
                with self.assertRaisesRegex(ValueError, 'modification time changed'):
                    frozen_input.stage_frozen_input(url, destination, expected, stage)
        with mock.patch.object(frozen_input.subprocess, 'run', return_value=SimpleNamespace(stdout='Size: 3\n')):
            with self.assertRaisesRegex(ValueError, 'lacks readable Size/MTime'):
                frozen_input.check_source_identity(url, expected)


if __name__ == '__main__':
    unittest.main()
