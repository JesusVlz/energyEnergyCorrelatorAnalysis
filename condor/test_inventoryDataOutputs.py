#!/usr/bin/env python3
"""Scoped metadata-inventory checks; all fixtures stay in temporary directories."""

from contextlib import redirect_stdout
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name('inventoryDataOutputs.py')
MODULE_SPEC = importlib.util.spec_from_file_location('inventoryDataOutputs', MODULE_PATH)
inventory = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(inventory)


def spec(sample, part, chunk=1, jobs=2):
    name = f'eecAnalysis_{sample}_mixedPerp_part{part:02d}_chunk{chunk:02d}_Run3_v2'
    return {'sample': sample, 'part': part, 'chunk': chunk,
            'request_name': name, 'expected_jobs': jobs, 'jobs': jobs}


class ScopeTests(unittest.TestCase):
    def setUp(self):
        self.specs = [spec('OO', 1), spec('OO', 2, 2, 3),
                      spec('ppRef', 1, jobs=4), spec('ppRef', 4, jobs=1)]

    def test_default_and_explicit_maximum_keep_all_tasks(self):
        for cutoff in (None, 2):
            with self.subTest(cutoff=cutoff):
                selected, scope = inventory.select_task_scope(self.specs, cutoff)
                self.assertEqual(selected, self.specs)
                self.assertTrue(scope['includes_all_authoritative_tasks'])
                self.assertEqual(scope['excluded_authoritative_request_names'], [])

    def test_cutoff_keeps_all_ppref_and_records_omitted_oo(self):
        selected, scope = inventory.select_task_scope(self.specs, 1)
        self.assertEqual(selected, [self.specs[0], *self.specs[2:]])
        self.assertEqual(scope['max_oo_part'], 1)
        self.assertEqual(scope['authoritative_max_oo_part'], 2)
        self.assertFalse(scope['includes_all_authoritative_tasks'])
        self.assertEqual(scope['selected_task_count'], 3)
        self.assertEqual(scope['authoritative_task_count'], 4)
        self.assertEqual(scope['excluded_authoritative_request_names'],
                         [self.specs[1]['request_name']])

    def test_invalid_cutoff_uses_authoritative_maximum(self):
        for cutoff in (-1, 0, 3):
            with self.subTest(cutoff=cutoff):
                with self.assertRaisesRegex(ValueError, 'between 1 and 2'):
                    inventory.select_task_scope(self.specs, cutoff)

    def test_inventory_preserves_full_counts_and_requires_same_recheck_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'source'
            source.mkdir()
            for sample, directory in [('OO', 'OOForestfiles'), ('ppRef', 'ppRefForestfiles')]:
                path = root / 'crab' / directory / 'manifest.json'
                path.parent.mkdir(parents=True)
                path.write_text(json.dumps({'units_per_job': 1,
                                            'tasks': [task for task in self.specs
                                                      if task['sample'] == sample]}))
            selected, unused = inventory.select_task_scope(self.specs, 1)
            for task in selected:
                directory = source / task['request_name'] / 'dataset' / 'submission' / '0000'
                directory.mkdir(parents=True)
                file = directory / (task['request_name'] + '_1.root')
                file.write_bytes(b'metadata-only fixture')
                os.utime(file, (time.time() - 600, time.time() - 600))
            output = root / 'inventory.json'
            arguments = [str(MODULE_PATH), '--repo-root', str(root), '--source-root', str(source),
                         '--output', str(output), '--max-oo-part', '1']
            with patch('sys.argv', arguments), redirect_stdout(io.StringIO()):
                inventory.main()
            report = json.loads(output.read_text())
            self.assertEqual(report['full_campaign_expected_jobs'], {'OO': 5, 'ppRef': 5})
            self.assertEqual(report['summary']['OO']['expected_jobs'], 2)
            self.assertEqual(report['summary']['OO']['full_campaign_expected_jobs'], 5)
            self.assertEqual(report['summary']['OO']['excluded_authoritative_jobs'], 3)
            self.assertEqual(report['summary']['ppRef']['expected_jobs'], 5)
            self.assertEqual(report['summary']['ppRef']['excluded_authoritative_jobs'], 0)
            excluded = report['selection_scope']['excluded_authoritative_request_names']
            self.assertEqual(excluded,
                             ['eecAnalysis_OO_mixedPerp_part02_chunk02_Run3_v3'])
            rechecked = root / 'rechecked.json'
            recheck_arguments = arguments[:]
            recheck_arguments[recheck_arguments.index('--output') + 1] = str(rechecked)
            recheck_arguments += ['--input-inventory', str(output)]
            with patch('sys.argv', recheck_arguments), redirect_stdout(io.StringIO()):
                inventory.main()
            self.assertEqual(json.loads(rechecked.read_text())['summary'], report['summary'])
            recheck_arguments = recheck_arguments[:-4] + recheck_arguments[-2:]
            recheck_arguments[recheck_arguments.index('--output') + 1] = str(root / 'mismatch.json')
            with patch('sys.argv', recheck_arguments), self.assertRaisesRegex(
                    ValueError, 'same --max-oo-part'):
                inventory.main()


class CandidateTests(unittest.TestCase):
    def setUp(self):
        self.source = Path('/example')
        self.task = spec('OO', 1)
        self.root = self.source / self.task['request_name']

    def record(self, job=1, namespace='submission', directory='0000'):
        name = self.task['request_name'] + f'_{job}.root'
        return {'path': str(self.root / 'dataset' / namespace / directory / name),
                'job_id': job, 'retry': None}

    def test_distinct_jobs_in_one_namespace_are_accepted(self):
        inventory.validate_candidates(self.task, [self.record(1), self.record(2)], self.source)

    def test_duplicate_job_across_numeric_directories_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate output job'):
            inventory.validate_candidates(self.task,
                                          [self.record(), self.record(directory='0001')],
                                          self.source)

    def test_retry_metadata_and_retry_filename_are_rejected(self):
        record = self.record()
        record['retry'] = 1
        with self.assertRaisesRegex(ValueError, 'Explicit retry'):
            inventory.validate_candidates(self.task, [record], self.source)
        record = self.record()
        record['path'] = record['path'].replace('_1.root', '_1_1.root')
        with self.assertRaisesRegex(ValueError, 'Unexpected output or retry filename'):
            inventory.validate_candidates(self.task, [record], self.source)

    def test_mixed_submission_namespaces_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Multiple submission namespaces'):
            inventory.validate_candidates(self.task,
                                          [self.record(1), self.record(2, namespace='other')],
                                          self.source)


if __name__ == '__main__':
    unittest.main()
