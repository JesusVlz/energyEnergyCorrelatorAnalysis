"""Verify combined-worker sequencing, complete coverage, provenance and failures."""
import contextlib
import gzip
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

spec = importlib.util.spec_from_file_location('combined', Path(__file__).with_name('pair_response_worker.py'))
worker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(worker)


class CombinedWorkerTest(unittest.TestCase):
    def fixture(self, root):
        analysis = root/'analysis'
        for name in ('eecAnalysis', 'validation/run3Pairs/batch_worker.py',
            'validation/run3Pairs/pair_counts.py', 'validation/run3Pairs/batch_pairs.py',
            'validation/run3Pairs/export_forest.py', 'validation/run3Unfolding/batch/response_worker.py',
            'validation/run3Unfolding/batch/cards/pp_R04_response.input'):
            path = analysis/name;path.parent.mkdir(parents=True, exist_ok=True);path.write_text(name)
        source = root/'staged.root';source.write_bytes(b'staged forest fixture')
        task = {'system': 'pp', 'file_index': 714, 'sample_namespace': 'physical-sample',
            'input': '/original/forest.root', 'entry_limit': 4,
            'input_identity': {'size_bytes': source.stat().st_size, 'entries': 4, 'sha256': worker.digest(source)}}
        manifest = root/'manifest.json';worker.save(manifest, {'tasks': [task]})
        args = SimpleNamespace(analysis_dir=analysis, manifest=manifest, input_local=source,
            output_dir=root/'products', scratch=root/'scratch', task_index=0)
        return args, task

    def helpers(self, args, task, failure=None):
        calls = []
        def invoke(command, cwd, check):
            calls.append(command)
            self.assertEqual(cwd, args.analysis_dir)
            self.assertTrue(check)
            if command[1].endswith('/batch_worker.py'):
                self.assertEqual(command[command.index('--input-local')+1], str(args.input_local))
                out = Path(command[command.index('--output-dir')+1]);out.mkdir()
                prefix = 'pp_714';numbers = {'raw': 4, 'selected': 3, 'train': 2, 'holdout': 1}
                counts = {'task': task, 'event_counts': numbers, 'raw_identity_scope': 'ALL_SCANNED_EVENTS',
                    'physics_core_sha256': worker.digest(args.analysis_dir/'validation/run3Pairs/pair_counts.py'),
                    'native_metadata': {'raw_input_sha256': worker.digest(args.input_local),
                        'frozen_task_manifest_sha256': worker.digest(args.manifest),
                        'native_binary_sha256': 'wrong' if failure == 'fingerprint' else worker.digest(args.analysis_dir/'eecAnalysis')}}
                with gzip.open(out/(prefix+'.counts.json.gz'), 'wt') as stream:json.dump(counts, stream)
                with gzip.open(out/(prefix+'.identities.bin.gz'), 'wb') as stream:
                    for index, marker in enumerate((1, 2, 0, 1)):
                        stream.write(worker.IDENTITY.pack(1, 2, index, marker))
                with gzip.open(out/(prefix+'.holdout.jsonl.gz'), 'wt') as stream:stream.write('{}\n')
                files = {}
                for kind, suffix in [('counts', 'counts.json.gz'), ('identities', 'identities.bin.gz'), ('holdout', 'holdout.jsonl.gz')]:
                    path = out/(prefix+'.'+suffix)
                    files[kind] = {'name': path.name, 'bytes': path.stat().st_size, 'sha256': worker.digest(path)}
                worker.save(out/(prefix+'.receipt.json'), {'task': task, 'event_counts': numbers, 'files': files})
            else:
                self.assertTrue(command[1].endswith('/response_worker.py'))
                self.assertEqual(command[2], str(args.input_local))
                self.assertEqual(command[command.index('--source-input')+1], task['input'])
                out = Path(command[3]);out.mkdir()
                (out/'train.log').write_text('retained native failure reason\n')
                if failure == 'native':raise subprocess.CalledProcessError(-6, command)
                ledger = out/'event_identity.csv.gz'
                with gzip.open(ledger, 'wt') as stream:
                    stream.write('source_entry,run,lumi,event,event_key,fold,fold_entry\n')
                    for index, fold in enumerate((0, 1, 1, 0)):
                        event = 123 if failure == 'identity' and index == 2 else index
                        stream.write('{0},1,2,{1},unused,{2},0\n'.format(index, event, fold))
                folds = {}
                for label in ('train', 'validation'):
                    path = out/(label+'_raw.root');path.write_bytes(label.encode())
                    folds[label] = {'bytes': path.stat().st_size, 'sha256': worker.digest(path)}
                receipt = {'status': 'COMPLETE_ENGINEERING_ONLY', 'skim': {
                    'input_entries': 4, 'processed_entries': 3 if failure == 'coverage' else 4,
                    'fold_entries': [2, 2], 'sample_namespace': task['sample_namespace'], 'ledger_sha256': worker.digest(ledger)},
                    'native_sha256': worker.digest(args.analysis_dir/'eecAnalysis'),
                    'pair_counts_sha256': worker.digest(args.analysis_dir/'validation/run3Pairs/pair_counts.py'),
                    'helper_sha256': worker.digest(args.analysis_dir/'validation/run3Unfolding/batch/response_worker.py'),
                    'card_sha256': worker.digest(args.analysis_dir/'validation/run3Unfolding/batch/cards/pp_R04_response.input'), 'folds': folds}
                worker.save(out/'receipt.json', receipt)
        return invoke, calls

    def test_success_binds_both_stages_to_same_input_and_writes_checked_markers_last(self):
        with tempfile.TemporaryDirectory() as temporary:
            args, task = self.fixture(Path(temporary));invoke, calls = self.helpers(args, task)
            with mock.patch.object(worker.subprocess, 'run', invoke):result = worker.run(args)
            self.assertEqual(len(calls), 2)
            self.assertEqual(result['shared_identity_audit']['matched_raw_physical_ids'], 4)
            receipt = worker.load(args.output_dir/'response/receipt.json')
            self.assertEqual(receipt['frozen_task'], task)
            self.assertEqual(receipt['campaign_manifest_sha256'], worker.digest(args.manifest))
            for stage in ('pairs', 'response'):
                directory = args.output_dir/stage;marker = worker.load(directory/'stage_complete.json')
                self.assertEqual(marker['status'], 'COMPLETE')
                self.assertEqual(set(marker['files']), {str(p.relative_to(directory)) for p in directory.rglob('*') if p.is_file() and p.name != 'stage_complete.json'})
                for name, record in marker['files'].items():
                    self.assertEqual(record['sha256'], worker.digest(directory/name))
                    self.assertEqual(record['bytes'], (directory/name).stat().st_size)

    def test_remote_input_is_bound_to_task_and_checked_before_and_after_staging(self):
        with tempfile.TemporaryDirectory() as temporary:
            args, task = self.fixture(Path(temporary))
            task['input_identity'].pop('sha256')
            task['input_identity']['mtime_ns'] = 10_000_000_000
            worker.save(args.manifest, {'tasks': [task]})
            data = args.input_local.read_bytes()
            args.input_local.unlink()
            args.input_url = 'root://test//original/forest.root'
            invoke, calls = self.helpers(args, task)
            def execute(command, cwd=None, check=True):
                if command[0] == 'xrdcp':
                    self.assertFalse(args.input_local.exists())
                    args.input_local.write_bytes(data)
                else:
                    invoke(command, cwd, check)
            metadata = {'size_bytes': len(data), 'mtime_seconds': 10}
            with mock.patch('frozen_input.source_metadata', return_value=metadata) as stat, mock.patch.object(worker.subprocess, 'run', execute):
                result = worker.run(args)
            self.assertEqual(stat.call_count, 2)
            self.assertEqual(len(calls), 2)
            self.assertEqual(result['source_metadata'], {'before': metadata, 'after': metadata})
            self.assertEqual(result['source_input_sha256'], worker.digest(args.input_local))

    def test_same_size_changed_staged_input_or_unverifiable_identity_stops_before_helpers(self):
        for missing_checksum in (False, True):
            with self.subTest(missing_checksum=missing_checksum), tempfile.TemporaryDirectory() as temporary:
                args, task = self.fixture(Path(temporary))
                if missing_checksum:
                    task['input_identity'].pop('sha256')
                    worker.save(args.manifest, {'tasks': [task]})
                else:
                    data = args.input_local.read_bytes()
                    args.input_local.write_bytes(b'X' + data[1:])
                with mock.patch.object(worker.subprocess, 'run') as invoke:
                    with self.assertRaisesRegex(ValueError, 'frozen SHA256|checksum differs'):
                        worker.run(args)
                    invoke.assert_not_called()
                self.assertFalse(list(args.output_dir.rglob('stage_complete.json')))

    def test_failure_keeps_logs_and_never_marks_partial_outputs_complete(self):
        for failure in ('native', 'coverage', 'identity', 'fingerprint'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as temporary:
                args, task = self.fixture(Path(temporary));invoke, calls = self.helpers(args, task, failure)
                stderr = io.StringIO()
                with mock.patch.object(worker.subprocess, 'run', invoke), contextlib.redirect_stderr(stderr):
                    with self.assertRaises((ValueError, subprocess.CalledProcessError)):worker.run(args)
                self.assertFalse(list(args.output_dir.rglob('stage_complete.json')))
                self.assertEqual(worker.load(args.output_dir/'failure.json')['status'], 'FAILED')
                if failure == 'fingerprint':self.assertEqual(len(calls), 1)
                else:self.assertIn('retained native failure reason', stderr.getvalue())


if __name__ == '__main__':
    unittest.main()
