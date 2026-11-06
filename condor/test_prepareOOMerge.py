"""Contract and preservation tests for frozen OO merge campaigns (no ROOT required)."""
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import types
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent



# The command-line preparer imports its worker from its own directory.
sys.path.insert(0, str(HERE))
try:
    import merge_closure_worker as worker
    import prepareOOMerge as prepare
finally:
    sys.path.pop(0)


class MergeFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='oo_merge_contract_')
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.production = self.base / 'production'
        self.output = self.base / 'outputs'
        for folder in (self.production / 'cards', self.production / 'mixing_pool', self.output / 'GenGen'):
            folder.mkdir(parents=True)
        self.source_dir = self.output / 'GenGen'
        self.url = lambda path: 'root://fixture//' + str(path).lstrip('/')
        self.patch_url = patch.object(prepare, 'eos_url', side_effect=self.url)
        self.patch_url.start()
        self.addCleanup(self.patch_url.stop)
        self.write(self.production / 'cards/GenGen.input', b'BackgroundMethods 6\n')
        self.write(self.production / 'analysis_payload.tar.gz', b'frozen analysis archive')
        provenance = {'source_commit': 'fixture', 'executable_sha256': 'e' * 64}
        self.write_json(self.production / 'campaign_provenance.json', provenance)
        self.sources = [{'path': '/eos/fixture/source%d.root' % index,
                         'source_entries': index + 2, 'root_uuid': 'source-%d' % index,
                         'size_bytes': index * 100} for index in (1, 2)]
        self.write_json(self.production / 'mixing_pool/signal_headers.manifest.json', {'files': self.sources})
        self.write_json(self.production / 'mixing_pool/source_audit.manifest.json', {'status': 'fixture_audit'})
        pool = {'signal_headers_sha256': worker.sha256(self.production / 'mixing_pool/signal_headers.manifest.json'),
                'source_audit_sha256': worker.sha256(self.production / 'mixing_pool/source_audit.manifest.json')}
        self.write_json(self.production / 'mixing_pool/manifest.json', pool)
        self.contract = {'mode': 'GenGen', 'data_type': 10,
                         'archive_sha256': worker.sha256(self.production / 'analysis_payload.tar.gz'),
                         'frozen_card_sha256': worker.sha256(self.production / 'cards/GenGen.input'),
                         'pool_manifest_sha256': worker.sha256(self.production / 'mixing_pool/manifest.json'),
                         'campaign_provenance_sha256': worker.sha256(self.production / 'campaign_provenance.json'),
                         'executable_sha256': provenance['executable_sha256'],
                         'signal_headers_sha256': pool['signal_headers_sha256'],
                         'source_audit_sha256': pool['source_audit_sha256'],
                         'partition_policy': 1, 'seed': 20261030, 'pool_cap': 100000,
                         'holdout_policy': 'fixed_reserved_donors_exclude_overlapping_signals',
                         'campaign_provenance': provenance}
        rows = []
        for index, source in enumerate(self.sources, 1):
            stem = 'sample_%d' % index
            input_url = self.url(source['path'])
            output_url = self.url(self.source_dir / (stem + '.root'))
            artifacts = {}
            for suffix in worker.SUFFIXES:
                # Distinct pool reuse counts must survive as separate job ledgers.
                content = ('%s%s: original per-job bytes %d\n' % (stem, suffix, index)).encode()
                path = self.source_dir / (stem + suffix)
                self.write(path, content)
                artifacts[path.name] = {'sha256': worker.sha256(path), 'size_bytes': len(content)}
            record = dict(self.contract, status=worker.PASS, input_url=input_url, output_url=output_url,
                          source=source, source_identity={'original_path': source['path'],
                                                         'root_uuid': source['root_uuid'],
                                                         'full_event_entries': source['source_entries']},
                          root_validation={'all_events': source['source_entries'],
                                           'inclusive_jet_entries': index, 'eec_entries': index * 2},
                          ledger_validation={'source_entries': source['source_entries'], 'pool_entries': 100000,
                                             'selected_signal_events': index, 'matched_signal_events': index,
                                             'matched_jets': index, 'contribution_rows': index * 2,
                                             'reserved_signal_failures': 0, 'partner_assignments': index * 2},
                          commands=[{'argv': ['analysis'], 'exit_code': 0}], artifacts=artifacts,
                          published_artifacts=[output_url[:-5] + suffix for suffix in worker.SUFFIXES])
            self.write_json(self.source_dir / (stem + '.run.json'), record)
            rows.append('\t'.join(['GenGen', input_url, output_url, stem, '4000']))
        self.write(self.production / 'jobs.tsv', ('\n'.join(rows) + '\n').encode())
        self.write_json(self.production / 'manifest.json', {
            'run_config': {'sample': 'OOMC', 'modes': ['GenGen'], 'background_methods': 6,
                           'data_type': 10, 'output_dir': str(self.output), 'tag': 'fixture_integrated'},
            'provenance': {'archive_sha256': self.contract['archive_sha256']},
            'inputs': [dict(row, url=self.url(row['path'])) for row in self.sources],
            'input_count': 2, 'job_count': 2})

    @staticmethod
    def write(path, data):
        path.write_bytes(data)

    @staticmethod
    def write_json(path, data):
        path.write_text(json.dumps(data, indent=2) + '\n')

    def audit(self):
        return prepare.audit_production(self.production)

    def record_path(self, index=1):
        return self.source_dir / ('sample_%d.run.json' % index)

    def mutate_record(self, mutate):
        path = self.record_path()
        record = json.loads(path.read_text())
        mutate(record)
        self.write_json(path, record)


class ProductionAuditTests(MergeFixture):
    def test_complete_production_binds_each_input_and_joint_exposure(self):
        manifest = self.audit()
        self.assertEqual(manifest['input_count'], 2)
        self.assertEqual(manifest['expected_root_totals'], {
            'all_events': 7, 'inclusive_jet_entries': 3, 'eec_entries': 6})
        self.assertEqual(manifest['ledger_totals']['partner_assignments'], 6)
        self.assertEqual(manifest['contract'], self.contract)
        self.assertEqual([entry['stem'] for entry in manifest['entries']], ['sample_1', 'sample_2'])
        for entry in manifest['entries']:
            self.assertEqual(entry['run_metadata']['sha256'], worker.sha256(self.source_dir / (entry['stem'] + '.run.json')))

    def test_incomplete_queue_unexpected_output_and_missing_sidecar_fail(self):
        queue = self.production / 'jobs.tsv'
        original = queue.read_bytes()
        queue.write_text(queue.read_text().splitlines()[0] + '\n')
        with self.assertRaisesRegex(ValueError, 'Incomplete production queue'):
            self.audit()
        queue.write_bytes(original)
        extra = self.source_dir / 'unexpected.root'
        extra.write_bytes(b'extra')
        with self.assertRaisesRegex(ValueError, 'Unexpected/missing ROOT output'):
            self.audit()
        extra.unlink()
        (self.source_dir / 'sample_1.pool.csv').unlink()
        with self.assertRaises(FileNotFoundError):
            self.audit()

    def test_mutated_frozen_payload_and_failed_job_are_rejected(self):
        archive = self.production / 'analysis_payload.tar.gz'
        original = archive.read_bytes()
        archive.write_bytes(b'changed archive')
        with self.assertRaisesRegex(ValueError, 'Frozen payload checksum mismatch'):
            self.audit()
        archive.write_bytes(original)
        self.mutate_record(lambda record: record.update(status='FAILED'))
        with self.assertRaisesRegex(ValueError, 'Unsuccessful run'):
            self.audit()

    def test_prepare_copies_helpers_and_binds_manifest_before_submission(self):
        workspace = self.base / 'merge_campaign'
        destination = self.base / 'merged'
        manifest = prepare.prepare(self.production, workspace, destination)
        frozen = json.loads((workspace / 'manifest.json').read_text())
        self.assertEqual(frozen, manifest)
        for name, digest in manifest['helper_sha256'].items():
            self.assertEqual(worker.sha256(workspace / name), digest)
        wrapper = (workspace / 'run_merge.sh').read_text()
        self.assertIn(worker.sha256(workspace / 'manifest.json'), wrapper)
        queue = (workspace / 'jobs.sub').read_text()
        self.assertIn('queue 1', queue)
        self.assertIn('transfer_output_files = ""', queue)
        self.assertEqual(list(destination.iterdir()), [])
        with self.assertRaisesRegex(ValueError, 'already exists'):
            prepare.prepare(self.production, workspace, destination)


class ContractValidationTests(MergeFixture):
    def test_record_contract_identity_exposure_and_publication_are_bound(self):
        manifest = self.audit()
        entry = manifest['entries'][0]
        mutations = {
            'failed worker': lambda r: r['commands'][0].update(exit_code=1),
            'other campaign': lambda r: r.update(archive_sha256='0' * 64),
            'other frozen card': lambda r: r.update(frozen_card_sha256='1' * 64),
            'other donor pool': lambda r: r.update(pool_manifest_sha256='2' * 64),
            'other source UUID': lambda r: r['source'].update(root_uuid='different'),
            'other source identity UUID': lambda r: r['source_identity'].update(root_uuid='different'),
            'truncated signal': lambda r: r['source_identity'].update(full_event_entries=1),
            'ROOT exposure': lambda r: r['root_validation'].update(all_events=1),
            'jet exposure': lambda r: r['root_validation'].update(inclusive_jet_entries=999),
            'donor scan': lambda r: r['ledger_validation'].update(pool_entries=1),
            'missing ledger': lambda r: r['artifacts'].pop('sample_1.pool.csv'),
            'missing publication': lambda r: r['published_artifacts'].pop(),
            'invalid hash': lambda r: r['artifacts']['sample_1.root'].update(sha256='z' * 64),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                record = copy.deepcopy(entry['run_record'])
                mutate(record)
                with self.assertRaises(ValueError):
                    worker.validate_record(record, entry, manifest['contract'])

    def test_manifest_rejects_incomplete_duplicate_mixed_or_unsupported_inputs(self):
        original = self.audit()
        mutations = {
            'incomplete': lambda m: m['entries'].pop(),
            'duplicate ROOT path': lambda m: m['entries'][1].update(output_url=m['entries'][0]['output_url']),
            'duplicate source': lambda m: m['entries'][1].update(input_url=m['entries'][0]['input_url']),
            'mixed campaigns': lambda m: m['entries'][1].update(output_url='root://fixture//another/campaign/sample_2.root'),
            'unsupported mode': lambda m: m['contract'].update(mode='RecoReco'),
            'incorrect expected exposure': lambda m: m['expected_root_totals'].update(all_events=999),
            'unsafe stem': lambda m: m['entries'][0].update(stem='../sample_1'),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                manifest = copy.deepcopy(original)
                mutate(manifest)
                with self.assertRaises(ValueError):
                    worker.validate_manifest(manifest)

    def test_staged_bytes_are_checked_even_when_size_is_unchanged(self):
        path = self.source_dir / 'sample_1.root'
        expected = {'size_bytes': path.stat().st_size, 'sha256': worker.sha256(path)}
        worker.checked_file(path, expected)
        content = path.read_bytes()
        path.write_bytes(bytes([content[0] ^ 1]) + content[1:])
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            worker.checked_file(path, expected)
        path.write_bytes(content + b'extra')
        with self.assertRaisesRegex(ValueError, 'Size mismatch'):
            worker.checked_file(path, expected)


class WorkerPreservationTests(MergeFixture):
    def execute_worker(self, mutate_staged=None, publish=False):
        manifest = self.audit()
        manifest.update(destination_url='root://fixture//merged')
        manifest_path = self.base / 'merge_manifest.json'
        self.write_json(manifest_path, manifest)
        self.publications = []
        self.merge_calls = []

        def transfer(argv, check=False):
            self.assertEqual(argv[0], 'xrdcp')
            source, destination = argv[-2:]
            if source.startswith('root://'):
                local_source = Path('/' + source.split('//', 2)[-1].lstrip('/'))
                shutil.copy2(local_source, destination)
                if mutate_staged and str(destination).endswith(mutate_staged):
                    data = Path(destination).read_bytes()
                    Path(destination).write_bytes(bytes([data[0] ^ 1]) + data[1:])
            else:
                self.publications.append(destination)
            return subprocess.CompletedProcess(argv, 0)

        def merge(inputs, output, totals):
            self.merge_calls.append((inputs, output, totals))
            self.assertEqual(len(inputs), 2)
            self.assertEqual([path.read_bytes() for path in inputs],
                             [(self.source_dir / ('sample_%d.root' % i)).read_bytes() for i in (1, 2)])
            output.write_bytes(b'merged ROOT fixture')
            return totals

        helper = types.SimpleNamespace(merge_root=merge)
        with patch.dict(sys.modules, {'merge_root': helper}), patch.object(worker.subprocess, 'run', side_effect=transfer):
            return worker.run(manifest_path, self.base / 'work', publish=publish)

    def test_joint_ledgers_and_original_run_records_preserve_exact_per_job_bytes(self):
        report = self.execute_worker()
        self.assertEqual(report['status'], 'PASS_MERGE_SOFTWARE_AND_EXPOSURE_CHECKS')
        self.assertFalse(report['published'])
        self.assertFalse(self.publications)
        bundle = self.base / 'work/merged_OOMC_GenGen.ledgers.tar.gz'
        with tarfile.open(bundle) as archive:
            names = archive.getnames()
            self.assertEqual(len(names), 11)  # frozen manifest + 2 * (run record + four ledgers)
            self.assertIn('merge_manifest.json', names)
            for index in (1, 2):
                for suffix in ('.run.json', *worker.SUFFIXES[1:]):
                    name = 'sample_%d%s' % (index, suffix)
                    self.assertEqual(archive.extractfile(name).read(), (self.source_dir / name).read_bytes())
            self.assertFalse(any(name.endswith('.root') for name in names))
        for name, expected in report['artifacts'].items():
            worker.checked_file(self.base / 'work' / name, expected)

    def test_changed_staged_ledger_stops_before_merge_or_publication(self):
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            self.execute_worker(mutate_staged='.pool.csv', publish=True)
        self.assertFalse(self.merge_calls)
        self.assertFalse(self.publications)


if __name__ == '__main__':
    unittest.main()
