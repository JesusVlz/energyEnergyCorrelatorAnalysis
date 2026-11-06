"""Regression tests for reduction resource sizing from transferred reports."""
import copy
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import prepareDataMerge as prepare


class ReductionResourceTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='eec-data-resource-')
        self.addCleanup(temporary.cleanup)
        self.workspace = Path(temporary.name)
        for name in ('nodes', 'reports', 'resources'):
            (self.workspace / name).mkdir()
        self.node = {'node_id': 'OO_part01_chunk01', 'kind': 'part',
                     'sample': 'OO', 'data_type': 9,
                     'parents': ['OO_batch001', 'OO_batch002'], 'inputs': []}
        self.reports = []
        for parent_id, size in zip(self.node['parents'], (8_000_000_000, 12_000_000_000)):
            url = 'root://eoscms.cern.ch//eos/cms/store/fixture/' + parent_id
            entry = {'url': url + '.root', 'report_url': url + '.json'}
            self.node['inputs'].append(entry)
            self.reports.append({'node_id': parent_id, 'sample': 'OO', 'data_type': 9,
                                 'status': 'PASS', 'artifact': {'size_bytes': size},
                                 'publication': {'root_url': entry['url'],
                                                 'report_url': entry['report_url'],
                                                 'root_published': True,
                                                 'report_published': True}})
        self.node_path = self.workspace / 'nodes' / (self.node['node_id'] + '.json')
        self.resource_path = self.workspace / 'resources' / (self.node['node_id'] + '.sub')
        self.resource_path.write_text('existing resource settings\n')
        self.write_fixture()

    def write_fixture(self):
        self.node_path.write_text(json.dumps(self.node))
        for parent_id, report in zip(self.node['parents'], self.reports):
            (self.workspace / 'reports' / (parent_id + '.json')).write_text(json.dumps(report))

    def run_update(self):
        with contextlib.redirect_stdout(io.StringIO()) as output:
            prepare.update_resources(self.node_path)
        return json.loads(output.getvalue())

    def assert_resources_unchanged(self):
        self.assertEqual(self.resource_path.read_text(), 'existing resource settings\n')
        self.assertFalse(self.resource_path.with_suffix('.tmp').exists())

    def test_sizes_part_and_sample_using_transferred_reports_only(self):
        # The remote fixture URLs do not exist; only campaign-local JSON is available.
        self.reports[1]['status'] = 'PARTIAL'
        for kind, memory in (('part', 16000), ('sample', 24000)):
            with self.subTest(kind=kind):
                self.node['kind'] = kind
                self.write_fixture()
                result = self.run_update()
                self.assertEqual(result['parent_bytes'], 20_000_000_000)
                self.assertEqual(result['resources']['disk_mb'], 80000)
                self.assertEqual(result['resources']['memory_mb'], memory)
                self.assertIn('request_disk = 80000MB\n', self.resource_path.read_text())
                self.assertFalse(self.resource_path.with_suffix('.tmp').exists())

    def test_missing_transferred_report_does_not_replace_resources(self):
        (self.workspace / 'reports' / (self.node['parents'][1] + '.json')).unlink()
        with self.assertRaises(FileNotFoundError):
            self.run_update()
        self.assert_resources_unchanged()

    def test_wrong_parent_identity_or_configuration_is_rejected(self):
        original = copy.deepcopy(self.reports[0])
        for key, value in (('node_id', 'OO_other_batch'), ('sample', 'ppRef'),
                           ('data_type', 11), ('status', 'FAIL')):
            with self.subTest(key=key):
                self.reports[0] = {**original, key: value}
                self.write_fixture()
                with self.assertRaises(ValueError):
                    self.run_update()
                self.assert_resources_unchanged()

    def test_wrong_publication_or_unpublished_parent_is_rejected(self):
        original = copy.deepcopy(self.reports[0])
        for key, value in (('root_url', 'root://wrong/input.root'),
                           ('report_url', 'root://wrong/input.json'),
                           ('root_published', False), ('report_published', False),
                           ('report_published', 'true')):
            with self.subTest(key=key, value=value):
                self.reports[0] = copy.deepcopy(original)
                self.reports[0]['publication'][key] = value
                self.write_fixture()
                with self.assertRaisesRegex(ValueError, 'Wrong/unpublished'):
                    self.run_update()
                self.assert_resources_unchanged()

    def test_invalid_artifact_sizes_are_rejected(self):
        for size in (0, -1, True, 1.5, '100'):
            with self.subTest(size=size):
                self.reports[0]['artifact']['size_bytes'] = size
                self.write_fixture()
                with self.assertRaisesRegex(ValueError, 'artifact size'):
                    self.run_update()
                self.assert_resources_unchanged()

    def test_parent_list_must_align_with_inputs(self):
        original = copy.deepcopy(self.node)
        for parents in ([], ['OO_batch001'], ['OO_batch001', 'OO_batch001'],
                        ['OO_batch002', 'OO_batch001'], ['../outside', 'OO_batch002']):
            with self.subTest(parents=parents):
                self.node = {**original, 'parents': parents}
                # Preserve the original reports while testing the frozen node mapping.
                self.node_path.write_text(json.dumps(self.node))
                with self.assertRaises(ValueError):
                    self.run_update()
                self.assert_resources_unchanged()


class SnapshotScopeTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='eec-data-scope-')
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.inventory = {'stable_snapshot': True, 'snapshot_started_utc': '2026-11-06T16:00:00+00:00',
                          'summary': {}, 'tasks': []}
        for sample in ('OO', 'ppRef'):
            directory = '/eos/cms/store/fixture/' + sample
            self.inventory['summary'][sample] = {'expected_jobs': 1, 'observed_files': 1}
            self.inventory['tasks'].append({'sample': sample, 'part': 1, 'chunk': 1,
                                           'request_name': sample + '_request', 'expected_jobs': 1,
                                           'missing_job_ids': [], 'directory': directory,
                                           'files': [{'path': directory + '/file_1.root',
                                                      'job_id': 1, 'bytes': 1000, 'mtime_ns': 123}]})

    def prepare_fixture(self):
        inventory_path = self.directory / 'inventory.json'
        inventory_path.write_text(json.dumps(self.inventory))
        workspace = self.directory / 'campaign'
        with mock.patch.object(prepare, 'eos_url', side_effect=lambda path: 'root://fixture/' + str(path)), \
                contextlib.redirect_stdout(io.StringIO()):
            summary = prepare.prepare(inventory_path, workspace, self.directory / 'output', 'fixture')
        submit = (workspace / 'submit.sh').read_text()
        self.assertIn('condor_submit_dag -no_submit -MaxJobs 20 merge.dag', submit)
        self.assertIn('condor_submit_dag -update_submit -MaxJobs 20 merge.dag', submit)
        snapshots = {sample: json.loads((workspace / 'nodes' / ('merged_' + sample + '.json')).read_text())['snapshot']
                     for sample in ('OO', 'ppRef')}
        return summary, snapshots, (workspace / 'README.txt').read_text()

    def test_complete_selected_subset_is_not_complete_authoritative_campaign(self):
        scope = {'max_oo_part': 1, 'authoritative_max_oo_part': 2,
                 'includes_all_authoritative_tasks': False,
                 'excluded_authoritative_request_names': ['OO_part02_request'],
                 'selected_task_count': 2, 'authoritative_task_count': 3}
        self.inventory['selection_scope'] = scope
        self.inventory['full_campaign_expected_jobs'] = {'OO': 2, 'ppRef': 1}
        summary, snapshots, readme = self.prepare_fixture()
        self.assertFalse(summary['snapshot_is_complete_campaign'])
        self.assertEqual(summary['selection_scope'], scope)
        self.assertEqual(summary['full_campaign_expected_jobs'], {'OO': 2, 'ppRef': 1})
        self.assertTrue(snapshots['OO']['complete'])
        self.assertEqual(snapshots['OO']['expected_jobs'], 1)
        self.assertEqual(snapshots['OO']['full_campaign_expected_jobs'], 2)
        self.assertEqual(snapshots['OO']['selection_scope'], scope)
        self.assertEqual(snapshots['ppRef']['selection_scope'], scope)
        self.assertIn('Requested OO coverage: parts 1-1 of 2 authoritative parts', readme)

    def test_older_inventory_defaults_full_counts_to_sample_summary(self):
        summary, snapshots, _ = self.prepare_fixture()
        self.assertTrue(summary['snapshot_is_complete_campaign'])
        self.assertNotIn('selection_scope', summary)
        self.assertEqual(summary['full_campaign_expected_jobs'], {'OO': 1, 'ppRef': 1})
        for snapshot in snapshots.values():
            self.assertEqual(snapshot['full_campaign_expected_jobs'], 1)
            self.assertNotIn('selection_scope', snapshot)


if __name__ == '__main__':
    unittest.main()
