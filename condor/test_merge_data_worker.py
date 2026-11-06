#!/usr/bin/env python3
"""Local ROOT regressions for CRAB data batching and parent reductions."""
from array import array
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import ROOT
import merge_data_worker as worker


class DataMergeWorkerTest(unittest.TestCase):
    def setUp(self):
        ROOT.gROOT.SetBatch(True)
        self.temporary = tempfile.TemporaryDirectory(prefix='eec-data-merge-')
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.sequence = 0

    def input(self, name, count=1, data_type=9, sumw2=True):
        path = self.directory / name
        output = ROOT.TFile(str(path), 'RECREATE')
        events = ROOT.TH1F('nEvents', 'Event selection', 9, -.5, 8.5)
        if sumw2:
            events.Sumw2()
        events.GetXaxis().SetTitle('selection')
        events.GetXaxis().SetBinLabel(1, 'all events')
        for index in range(1, 10):
            events.SetBinContent(index, count)
            if sumw2:
                events.SetBinError(index, count ** .5)
        events.SetEntries(count)
        events.Write()
        for name in ('inclusiveJet', 'energyEnergyCorrelator'):
            histogram = ROOT.TH1D(name, name, 2, 0, 2)
            histogram.Sumw2()
            histogram.Fill(.5)
            histogram.Write()
        ROOT.TObjString(str(path)).Write('InputFiles')
        card = output.mkdir('JCard')
        card.cd()
        for key, value in (('DataType', data_type), ('McCorrelationType', 0), ('BackgroundMethods', 6)):
            vector = ROOT.TVectorF(1, 1)
            vector.GetMatrixArray()[0] = value
            vector.Write(key)
        ROOT.TObjString('fixture-git').Write('GitHash')
        output.Close()
        return path

    def entry(self, path):
        return {'url': str(path), 'size_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns}

    def manifest(self, inputs, kind='batch', sample='OO'):
        self.sequence += 1
        identifier = 'node_%d' % self.sequence
        destination = self.directory / 'published' / identifier
        manifest = {'node_id': identifier, 'sample': sample, 'kind': kind,
                    'data_type': 9 if sample == 'OO' else 11, 'inputs': inputs,
                    'output_url': str(destination.with_suffix('.root')),
                    'report_url': str(destination.with_suffix('.json')),
                    'snapshot': {'expected_jobs': 4, 'observed_jobs': 3, 'missing_job_ids': [4]}}
        path = self.directory / (identifier + '.node.json')
        path.write_text(json.dumps(manifest))
        return path, self.directory / (identifier + '.work'), manifest

    def parent_entry(self, manifest):
        return {'url': manifest['output_url'], 'report_url': manifest['report_url']}

    def test_exact_counts_above_float_integer_precision_and_original_unchanged(self):
        first = self.input('large.root', 2 ** 24)
        second = self.input('one.root', 1)
        before = worker.sha256(first)
        manifest, work, expected = self.manifest([self.entry(first), self.entry(second)])
        report = worker.run(manifest, work, local_test=True, publish=False)
        self.assertEqual(report['status'], 'PASS')
        self.assertEqual(report['root_totals']['all_events'], 2 ** 24 + 1)
        self.assertEqual(report['event_counter_normalization']['promoted_input_count'], 2)
        self.assertEqual(report['snapshot'], expected['snapshot'])
        self.assertEqual(worker.sha256(first), before)
        output = ROOT.TFile.Open(str(work / 'result.root'))
        events = output.Get('nEvents')
        self.assertEqual(events.ClassName(), 'TH1D')
        self.assertEqual(events.GetBinContent(1), 2 ** 24 + 1)
        self.assertEqual(events.GetXaxis().GetBinLabel(1), 'all events')
        self.assertEqual(events.GetXaxis().GetTitle(), 'selection')
        self.assertEqual(events.GetEntries(), 2 ** 24 + 1)
        self.assertGreater(events.GetSumw2N(), 0)
        output.Close()

    def test_promotion_preserves_absent_sumw2(self):
        path = self.input('without_errors.root', 3, sumw2=False)
        self.assertTrue(worker._promote_event_counter(path))
        self.assertFalse(worker._promote_event_counter(path))
        file = ROOT.TFile.Open(str(path))
        self.assertEqual(file.Get('nEvents').GetSumw2N(), 0)
        file.Close()

    def test_zero_and_corrupt_raw_files_are_explicit_exclusions(self):
        good = self.input('good.root')
        empty = self.directory / 'empty.root'
        empty.touch()
        corrupt = self.directory / 'corrupt.root'
        corrupt.write_bytes(b'not a ROOT file')
        manifest, work, _ = self.manifest([self.entry(good), self.entry(empty), self.entry(corrupt)])
        report = worker.run(manifest, work, local_test=True, publish=False)
        self.assertEqual(report['status'], 'PARTIAL')
        self.assertEqual(report['included_source_urls'], [str(good)])
        self.assertEqual({row['url'] for row in report['excluded_files']}, {str(empty), str(corrupt)})
        self.assertEqual((work / 'good_files.txt').read_text().splitlines(), [str(good)])
        self.assertEqual(set((work / 'bad_files.txt').read_text().splitlines()), {str(empty), str(corrupt)})

    def test_all_bad_has_failure_report_and_diagnostic_artifacts(self):
        empty = self.directory / 'empty.root'
        empty.touch()
        manifest, work, _ = self.manifest([self.entry(empty)])
        with self.assertRaisesRegex(ValueError, 'No valid ROOT'):
            worker.run(manifest, work, local_test=True, publish=False)
        report = json.loads((work / 'result.json').read_text())
        self.assertEqual(report['status'], 'FAIL')
        self.assertIsNone(report['artifact'])
        self.assertEqual((work / 'result.root').stat().st_size, 0)
        self.assertTrue((work / 'good_files.txt').is_file())
        self.assertTrue((work / 'bad_files.txt').is_file())

    def test_wrong_card_and_changed_size_fail_instead_of_skipping(self):
        wrong = self.input('wrong_type.root', data_type=11)
        manifest, work, _ = self.manifest([self.entry(wrong)])
        with self.assertRaisesRegex(ValueError, 'DataType'):
            worker.run(manifest, work, local_test=True, publish=False)
        report = json.loads((work / 'result.json').read_text())
        self.assertEqual(report['status'], 'FAIL')
        self.assertEqual(report['excluded_files'], [])
        good = self.input('good.root')
        changed = self.entry(good)
        changed['size_bytes'] += 1
        manifest, work, _ = self.manifest([changed])
        with self.assertRaisesRegex(ValueError, 'size changed'):
            worker.run(manifest, work, local_test=True, publish=False)
        self.assertEqual(json.loads((work / 'result.json').read_text())['excluded_files'], [])

    def test_same_size_source_replacement_before_or_during_staging_fails(self):
        for during_staging in (False, True):
            with self.subTest(during_staging=during_staging):
                good = self.input('changed_%s.root' % during_staging)
                entry = self.entry(good)
                manifest, work, _ = self.manifest([entry])
                def replace_source():
                    os.utime(good, ns=(entry['mtime_ns'], entry['mtime_ns'] + 2_000_000_000))
                original_stage = worker._stage
                def stage(*args):
                    original_stage(*args)
                    replace_source()
                if not during_staging:
                    replace_source()
                with mock.patch.object(worker, '_stage', side_effect=stage), mock.patch.object(worker, '_check_data_card') as inspect:
                    with self.assertRaisesRegex(ValueError, 'modification time changed'):
                        worker.run(manifest, work, local_test=True, publish=False)
                    inspect.assert_not_called()
                self.assertEqual(good.stat().st_size, entry['size_bytes'])
                report = json.loads((work / 'result.json').read_text())
                self.assertEqual(report['status'], 'FAIL')
                self.assertEqual(report['excluded_files'], [])

    def test_missing_source_transport_fails_instead_of_skipping(self):
        manifest, work, _ = self.manifest([{'url': str(self.directory / 'missing.root'), 'size_bytes': 42, 'mtime_ns': 1}])
        with self.assertRaises(FileNotFoundError):
            worker.run(manifest, work, local_test=True, publish=False)
        report = json.loads((work / 'result.json').read_text())
        self.assertEqual(report['status'], 'FAIL')
        self.assertEqual(report['excluded_files'], [])

    def test_parent_reduction_propagates_exclusions_and_checks_overlap(self):
        good = self.input('good.root')
        empty = self.directory / 'empty.root'
        empty.touch()
        source, work, batch = self.manifest([self.entry(good), self.entry(empty)])
        worker.run(source, work, local_test=True)
        source, work, part = self.manifest([self.parent_entry(batch)], kind='part')
        report = worker.run(source, work, local_test=True)
        self.assertEqual(report['status'], 'PARTIAL')
        self.assertEqual(report['included_source_urls'], [str(good)])
        self.assertEqual([row['url'] for row in report['excluded_files']], [str(empty)])
        self.assertEqual(report['event_counter_normalization']['promoted_input_count'], 0)
        # The batch and its reduction represent the same original source files.
        source, work, _ = self.manifest([self.parent_entry(batch), self.parent_entry(part)], kind='sample')
        with self.assertRaisesRegex(ValueError, 'overlap'):
            worker.run(source, work, local_test=True)
        self.assertEqual(json.loads((work / 'result.json').read_text())['status'], 'FAIL')

    def test_disjoint_batches_reduce_through_multi_parent_parts_and_sample(self):
        originals = []
        batches = []
        for index in range(1, 4):
            original = self.input('source_%d.root' % index, count=index)
            originals.append(str(original))
            source, work, batch = self.manifest([self.entry(original)])
            worker.run(source, work, local_test=True)
            batches.append(batch)
        parts = []
        for group in (batches[:2], batches[2:]):
            source, work, part = self.manifest([self.parent_entry(batch) for batch in group], kind='part')
            report = worker.run(source, work, local_test=True)
            self.assertEqual(report['status'], 'PASS')
            parts.append(part)
        source, work, sample = self.manifest([self.parent_entry(part) for part in parts], kind='sample')
        report = worker.run(source, work, local_test=True)
        self.assertEqual(report['status'], 'PASS')
        self.assertEqual(report['included_source_urls'], originals)
        self.assertEqual(report['included_source_count'], 3)
        self.assertEqual(report['root_totals']['all_events'], 6)
        self.assertEqual(report['root_totals']['inclusive_jet_entries'], 3)
        self.assertEqual(report['root_totals']['eec_entries'], 3)
        self.assertEqual(report['event_counter_normalization']['promoted_input_count'], 0)
        final = worker.root_merge.inspect_root(sample['output_url'])
        self.assertEqual(final['n_events_bins'], [0] + [6] * 9 + [0])
        self.assertEqual(final['card']['DataType']['value'], [9])

    def test_parent_hash_is_checked_before_any_staged_modification(self):
        good = self.input('good.root')
        source, work, batch = self.manifest([self.entry(good)])
        worker.run(source, work, local_test=True)
        with Path(batch['output_url']).open('ab') as stream:
            stream.write(b'modified')
        source, work, _ = self.manifest([self.parent_entry(batch)], kind='part')
        with mock.patch.object(worker, '_promote_event_counter') as promote:
            with self.assertRaisesRegex(ValueError, 'size/SHA256'):
                worker.run(source, work, local_test=True)
            promote.assert_not_called()

    def test_existing_publication_is_never_overwritten(self):
        good = self.input('good.root')
        source, work, manifest = self.manifest([self.entry(good)])
        destination = Path(manifest['output_url'])
        destination.parent.mkdir()
        destination.write_bytes(b'existing result')
        with self.assertRaisesRegex(ValueError, 'already exists'):
            worker.run(source, work, local_test=True)
        self.assertEqual(destination.read_bytes(), b'existing result')
        self.assertFalse(Path(manifest['report_url']).exists())
        self.assertEqual(json.loads((work / 'result.json').read_text())['status'], 'FAIL')


if __name__ == '__main__':
    unittest.main()
