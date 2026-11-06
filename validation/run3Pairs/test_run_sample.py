"""A resumed bounded sample must retain the original request and actual bytes."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import run_sample


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.source = self.folder/'input.root'
        self.source.write_bytes(b'input fixture')
        self.output = self.folder/'events.jsonl'
        self.output.write_text('{"id": 1}\n')
        self.task = dict(system='OO', input=str(self.source), output=str(self.output),
                         command=['export_forest.py', '--max-events', '10'])
        self.sources = {'analyzer': 'source hash'}
        self.metadata = {'counts': {'exported_events': 1},
                         'output_sha256': run_sample.digest(self.output)}
        self.metadata_path = Path(str(self.output)+'.metadata.json')
        self.metadata_path.write_text(json.dumps(self.metadata))
        self.saved = dict(self.task, exit_code=0, source_sha256=self.sources,
                          input_sha256=run_sample.digest(self.source), metadata=self.metadata)
        self.output.with_suffix('.run.json').write_text(json.dumps(self.saved))

    def reusable(self, task=None, sources=None):
        return run_sample.reusable_result(task or self.task, sources or self.sources)

    def test_unchanged_export_is_reused(self):
        self.assertEqual(self.reusable(), self.saved)

    def test_changed_event_cap_or_source_is_not_reused(self):
        task = copy.deepcopy(self.task)
        task['command'][-1] = '20'
        self.assertIsNone(self.reusable(task=task))
        self.assertIsNone(self.reusable(sources={'analyzer': 'changed'}))

    def test_changed_input_is_not_reused(self):
        self.source.write_bytes(b'changed input')
        self.assertIsNone(self.reusable())

    def test_changed_or_missing_export_is_not_reused(self):
        self.output.write_text('{"id": 2}\n')
        self.assertIsNone(self.reusable())
        self.output.unlink()
        self.assertIsNone(self.reusable())

    def test_changed_or_missing_metadata_is_not_reused(self):
        self.metadata_path.write_text('{}')
        self.assertIsNone(self.reusable())
        self.metadata_path.unlink()
        self.assertIsNone(self.reusable())


if __name__ == '__main__':
    unittest.main()
