"""Campaign safety checks using tiny fixtures; never run native MC."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import all_mc_campaign as campaign


class CampaignTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = Path(self.temp.name)/'repo'
        self.repo.mkdir()
        self.inputs = Path(self.temp.name)/'inputs'
        self.inputs.mkdir()
        self.output = self.repo/'results'
        self.sources = []
        for name in ('eecAnalysis', 'src/EECAnalyzer.cxx', 'validation/run3Pairs/export_forest.py', 'calibration.root'):
            path = self.repo/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('frozen '+name+'\n')
            self.sources.append(path)
        self.payload = self.repo/'calibration.root'
        self.samples = {'OO': (str(self.inputs), 'forest_', 'existing-namespace')}
        self.headers = lambda path: int(path.read_text())

    def prepare(self, counts=((1, 2),)):
        for number, entries in counts:
            (self.inputs/f'forest_{number}.root').write_text(str(entries))
        with contextlib.redirect_stdout(io.StringIO()):
            return campaign.prepare(self.output, samples=self.samples, sources=self.sources,
                                    entry_reader=self.headers, repo=self.repo)

    def fake_export(self, command, **kwargs):
        source, target = Path(command[2]), Path(command[3])
        entries = int(command[-1])
        work = target.parent/(target.stem+'.native')
        work.mkdir()
        target.write_text('\n'.join(json.dumps({'event': i}) for i in range(entries))+'\n')
        for population in ('reco', 'gen'):
            (work/(population+'.jsonl.events.jsonl')).write_text(''.join(
                json.dumps({'source_entry': i, 'id': [1, 1, i]})+'\n' for i in range(entries)))
            for suffix in ('.root', '.input', '.jsonl', '.log'):
                (work/(population+suffix)).write_text('fixture\n')
        metadata = dict(system='OO', input=str(source), sample_namespace='existing-namespace',
                        input_size_bytes=source.stat().st_size, entry_limit_per_file=entries,
                        tool_sha256=campaign.digest(self.repo/'validation/run3Pairs/export_forest.py'),
                        analyzer_sha256=campaign.digest(self.repo/'src/EECAnalyzer.cxx'),
                        payload=str(self.payload), payload_sha256=campaign.digest(self.payload),
                        output_sha256=campaign.digest(target), native_outputs=str(work),
                        counts={'exported_events': entries, 'selected_jets_R04': 1})
        Path(str(target)+'.metadata.json').write_text(json.dumps(metadata))
        return SimpleNamespace(returncode=0)

    def execute(self, resume=False, effect=None):
        with patch.object(campaign, 'root_entries', self.headers), \
             patch.object(campaign.subprocess, 'run', side_effect=effect or self.fake_export) as process, \
             contextlib.redirect_stdout(io.StringIO()):
            result = campaign.execute(self.output, workers=1, resume=resume)
            return result, process.call_count

    def test_real_root_header_entry_count(self):
        import ROOT
        from array import array
        path = self.inputs/'tiny_header.root'
        source = ROOT.TFile(str(path), 'RECREATE')
        try:
            source.mkdir('hiEvtAnalyzer').cd()
            tree = ROOT.TTree('HiTree', 'Tiny header fixture, no MC')
            value = array('i', [0])
            tree.Branch('entry', value, 'entry/I')
            for value[0] in range(3):
                tree.Fill()
            tree.Write()
        finally:
            source.Close()
        self.assertEqual(campaign.root_entries(path), 3)

    def test_real_root_missing_event_tree_rejected(self):
        import ROOT
        path = self.inputs/'missing_tree.root'
        source = ROOT.TFile(str(path), 'RECREATE')
        source.Close()
        with self.assertRaisesRegex(ValueError, 'Missing hiEvtAnalyzer/HiTree'):
            campaign.root_entries(path)

    def test_discovery_uses_real_files_across_missing_index(self):
        self.prepare(((187, 2), (189, 3), (334, 1)))
        manifest, _ = campaign.load_manifest(self.output)
        self.assertEqual([t['file_index'] for t in manifest['tasks']], [187, 189, 334])
        self.assertIn(188, manifest['discovery']['OO']['missing_indices'])
        self.assertEqual(manifest['expected']['OO'], {'files': 3, 'raw_entries': 6})
        self.assertEqual([t['command'][-1] for t in manifest['tasks']], ['2', '3', '1'])
        self.assertIsNone(manifest['events_per_file'])
        self.assertEqual(len((self.output/'OO_inputs.txt').read_text().splitlines()), 3)

    def test_prepare_and_status_never_launch(self):
        with patch.object(campaign.subprocess, 'run') as process:
            self.prepare()
            with contextlib.redirect_stdout(io.StringIO()):
                summary = campaign.status(self.output)
            process.assert_not_called()
        self.assertEqual(summary['status'], 'INCOMPLETE')
        self.assertEqual(summary['pending_files'], 1)
        self.assertEqual(summary['systems']['OO']['scanned_entries'], 0)

    def test_incomplete_status_cli_blocks_derivation_without_launching(self):
        self.prepare()
        with patch.object(campaign.sys, 'argv', ['all_mc_campaign.py', 'status', str(self.output)]), \
             patch.object(campaign.subprocess, 'run') as process, \
             contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as stopped:
                campaign.main()
            self.assertEqual(stopped.exception.code, 2)
            process.assert_not_called()

    def test_resume_skips_only_validated_success(self):
        self.prepare(((1, 2), (2, 3)))
        first, calls = self.execute()
        self.assertEqual((first['status'], calls), ('COMPLETE', 2))
        self.assertEqual(first['systems']['OO']['scanned_entries'], 5)
        second, calls = self.execute(resume=True)
        self.assertEqual((second['status'], calls), ('COMPLETE', 0))

    def test_changed_input_prevents_any_launch(self):
        self.prepare()
        (self.inputs/'forest_1.root').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'Changed frozen'):
            self.execute()

    def test_missing_input_prevents_launch(self):
        self.prepare()
        (self.inputs/'forest_1.root').unlink()
        with self.assertRaises(FileNotFoundError):
            self.execute()

    def test_changed_source_prevents_launch(self):
        self.prepare()
        self.sources[0].write_text('changed binary')
        with self.assertRaisesRegex(ValueError, 'Changed frozen'):
            self.execute()

    def test_missing_successful_export_is_not_silently_rerun(self):
        self.prepare()
        self.execute()
        (self.output/'OO_001.jsonl').unlink()
        with self.assertRaises(FileNotFoundError):
            self.execute(resume=True)

    def test_changed_output_is_not_silently_rerun(self):
        self.prepare()
        self.execute()
        (self.output/'OO_001.jsonl').write_text('changed output')
        with self.assertRaisesRegex(ValueError, 'Export SHA256 mismatch'):
            self.execute(resume=True)

    def test_changed_metadata_rejected(self):
        self.prepare()
        self.execute()
        path = self.output/'OO_001.jsonl.metadata.json'
        metadata = json.loads(path.read_text())
        metadata['counts']['selected_jets_R04'] = 99
        path.write_text(json.dumps(metadata))
        with self.assertRaisesRegex(ValueError, 'export or metadata changed'):
            self.execute(resume=True)

    def test_failed_attempt_is_preserved_on_explicit_resume(self):
        self.prepare()
        failed, calls = self.execute(effect=lambda *a, **k: SimpleNamespace(returncode=13))
        self.assertEqual((len(failed['failed']), calls), (1, 1))
        with self.assertRaisesRegex(ValueError, 'requires --resume'):
            self.execute()
        successful, calls = self.execute(resume=True)
        self.assertEqual((successful['status'], calls), ('COMPLETE', 1))
        archives = list((self.output/'failures').rglob('OO_001.run.json'))
        self.assertEqual(len(archives), 1)
        self.assertEqual(json.loads(archives[0].read_text())['exit_code'], 13)

    def test_success_process_with_short_ledger_is_failure(self):
        self.prepare()
        def short(command, **kwargs):
            response = self.fake_export(command, **kwargs)
            target = Path(command[3])
            (target.parent/(target.stem+'.native')/'reco.jsonl.events.jsonl').write_text('')
            return response
        summary, _ = self.execute(effect=short)
        self.assertEqual(len(summary['failed']), 1)
        saved = json.loads((self.output/'OO_001.run.json').read_text())
        self.assertIn('full frozen input entries', saved['error'])

    def test_unreceipted_partial_export_requires_inspection(self):
        self.prepare()
        (self.output/'OO_001.jsonl').write_text('partial')
        with self.assertRaisesRegex(ValueError, 'Unreceipted partial'):
            self.execute(resume=True)


if __name__ == '__main__':
    unittest.main()
