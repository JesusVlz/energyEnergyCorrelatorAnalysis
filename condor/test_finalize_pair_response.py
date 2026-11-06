"""Archive integrity, physical-task binding and complete-coverage gate tests."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('finalizer', Path(__file__).with_name('finalize_pair_response.py'))
finalizer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(finalizer)


def task(index):
    return {'system': 'pp', 'file_index': index, 'input': '/forest/input_{}.root'.format(index),
            'entry_limit': 4, 'input_identity': {'entries': 4}, 'sample_namespace': 'fixture'}


def fixture(path, current_task, manifest_hash, corrupt=False, wrong_task=False, extra_member=None):
    provenance = {'status': 'COMPLETE', 'frozen_task': current_task,
                  'campaign_manifest_sha256': manifest_hash,
                  'shared_identity_audit': {'status': 'PASS', 'matched_raw_physical_ids': 4},
                  'source_fingerprints': {'eecAnalysis': 'f' * 64}, 'driver_sha256': 'e' * 64}
    pair_receipt = {'task': current_task, 'event_counts': {'raw': 4}}
    response_receipt = {'status': 'COMPLETE_ENGINEERING_ONLY', 'frozen_task': current_task,
                        'campaign_manifest_sha256': manifest_hash,
                        'skim': {'input_entries': 4, 'processed_entries': 4}}
    with tarfile.open(path, 'w:gz') as archive:
        for kind, receipt_name, receipt in [('pairs', finalizer.prefix(current_task) + '.receipt.json', pair_receipt),
                                          ('response', 'receipt.json', response_receipt)]:
            files = {receipt_name: json.dumps(receipt).encode(),
                     'combined_worker_provenance.json': json.dumps(provenance).encode(),
                     'payload.bin': b'valid physics artifact'}
            marker = {'status': 'COMPLETE', 'frozen_task': task(999) if wrong_task else current_task,
                      'campaign_manifest_sha256': manifest_hash,
                      'files': {name: {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
                                for name, data in files.items()}}
            files['stage_complete.json'] = json.dumps(marker).encode()
            if corrupt:
                files['payload.bin'] = b'wrong physics artifact'
            for name, data in files.items():
                member = tarfile.TarInfo(kind + '/' + name)
                member.size = len(data)
                archive.addfile(member, io.BytesIO(data))
        if extra_member:
            archive.addfile(extra_member)


class FinalizerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.archive = self.root / 'input.tar.gz'

    def tearDown(self):
        self.temporary.cleanup()

    def test_extracts_only_requested_kind_with_task_provenance(self):
        fixture(self.archive, task(1), 'a' * 64)
        destination = self.root / 'selected'
        provenance = finalizer.extract_checked(self.archive, destination, 'pairs', task(1), 'a' * 64)
        self.assertEqual(provenance['shared_identity_audit']['matched_raw_physical_ids'], 4)
        self.assertTrue((destination / 'pp_001.receipt.json').is_file())
        self.assertFalse((destination / 'receipt.json').exists())
        self.assertFalse((destination / 'response').exists())

    def test_changed_artifact_rejected(self):
        fixture(self.archive, task(1), 'a' * 64, corrupt=True)
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            finalizer.extract_checked(self.archive, self.root / 'selected', 'pairs', task(1), 'a' * 64)

    def test_wrong_frozen_task_or_manifest_rejected(self):
        fixture(self.archive, task(1), 'a' * 64, wrong_task=True)
        with self.assertRaisesRegex(ValueError, 'frozen archive task/manifest'):
            finalizer.extract_checked(self.archive, self.root / 'wrong_task', 'pairs', task(1), 'a' * 64)
        fixture(self.archive, task(1), 'a' * 64)
        with self.assertRaisesRegex(ValueError, 'frozen archive task/manifest'):
            finalizer.extract_checked(self.archive, self.root / 'wrong_manifest', 'pairs', task(1), 'b' * 64)

    def test_traversal_or_link_rejected_even_in_unselected_kind(self):
        for index, name in enumerate(('response/../../escape', 'response/link')):
            member = tarfile.TarInfo(name)
            if index:
                member.type = tarfile.SYMTYPE
                member.linkname = '../../escape'
            fixture(self.archive, task(1), 'a' * 64, extra_member=member)
            with self.assertRaises(ValueError):
                finalizer.extract_checked(self.archive, self.root / ('bad' + str(index)), 'pairs', task(1), 'a' * 64)

    def test_duplicate_member_rejected(self):
        fixture(self.archive, task(1), 'a' * 64, extra_member=tarfile.TarInfo('pairs/payload.bin'))
        with self.assertRaisesRegex(ValueError, 'Duplicate archive member'):
            finalizer.extract_checked(self.archive, self.root / 'bad', 'pairs', task(1), 'a' * 64)

    def prepare_run(self):
        manifest = self.root / 'manifest.json'
        finalizer.save(manifest, {'tasks': [task(1), task(2)], 'counts': {'pp': 2},
                                  'outputs': {'pp': 'root://example//outputs/pp'}})
        archives = {}
        for current_task in [task(1), task(2)]:
            name = finalizer.prefix(current_task)
            archives[name] = self.root / (name + '.tar.gz')
            fixture(archives[name], current_task, finalizer.digest(manifest))
        args = SimpleNamespace(manifest=manifest, system='pp', kind='response',
                               work_dir=self.root / 'work', output_dir=self.root / 'output')
        return args, archives

    def test_missing_expected_archive_prevents_any_computation(self):
        args, archives = self.prepare_run()
        def fetch(url, destination):
            if 'pp_002' in url:
                raise subprocess.CalledProcessError(54, ['xrdcp'])
            shutil.copyfile(archives['pp_001'], destination)
        with patch.object(finalizer, 'fetch', side_effect=fetch), patch.object(finalizer, 'run_pipeline') as pipeline:
            with self.assertRaises(subprocess.CalledProcessError):
                finalizer.finalize(args)
            pipeline.assert_not_called()
        self.assertFalse(args.output_dir.exists())

    def test_all_archives_verified_before_pipeline_and_receipt(self):
        args, archives = self.prepare_run()
        fetched = []
        def fetch(url, destination):
            name = Path(url).name.removesuffix('.tar.gz')
            fetched.append(name)
            shutil.copyfile(archives[name], destination)
        def pipeline(manifest_path, manifest, system, kind, tasks, folders, work, output):
            self.assertEqual(fetched, ['pp_001', 'pp_002'])
            self.assertEqual(finalizer.load(work / 'verified_inputs.json')['raw_entries'], 8)
            self.assertEqual(len(folders), 2)
            self.assertFalse(list(work.glob('*.tar.gz')))
            output.mkdir()
            (output / 'result.root').write_bytes(b'fixture')
            return [['existing-merge-helper']]
        with patch.object(finalizer, 'fetch', side_effect=fetch), patch.object(finalizer, 'run_pipeline', side_effect=pipeline):
            report = finalizer.finalize(args)
        self.assertEqual(report['files'], 2)
        self.assertEqual(report['raw_entries'], 8)
        self.assertEqual(report['outputs']['result.root']['sha256'], hashlib.sha256(b'fixture').hexdigest())
        self.assertTrue((args.output_dir / 'finalization_receipt.json').is_file())


if __name__ == '__main__':
    unittest.main()
