#!/usr/bin/env python3
"""Merge frozen CRAB data snapshots in scratch, with explicit corrupt-file exclusions.

Only zero-byte, unreadable, zombie, or recovered raw ROOT files may be excluded.
Transport, size, configuration, schema, histogram, and parent provenance failures
abort the node. Successful reports preserve original source membership so that
later sample merges can reject overlapping parts.
"""
import argparse
from array import array
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import posixpath
import re
import shutil
import subprocess
import sys
from urllib.parse import urlsplit

import merge_root as root_merge
from frozen_input import stage_frozen_input

SUCCESS_STATUSES = ('PASS', 'PARTIAL')
SAMPLE_TYPES = {'OO': 9, 'ppRef': 11}
TOTAL_KEYS = ('all_events', 'inclusive_jet_entries', 'eec_entries')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    checksum = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def _location(value, local_test=False):
    require(isinstance(value, str) and bool(value) and
            not any(char in value for char in ('\n', '\r', '\x00')), 'Invalid file location')
    parsed = urlsplit(value)
    if parsed.scheme:
        require(parsed.scheme == 'root' and bool(parsed.netloc) and
                parsed.path.startswith('/') and not parsed.query and not parsed.fragment,
                'Expected an absolute XRootD file URL: ' + value)
        return ('root://' + parsed.netloc, '/' + parsed.path.lstrip('/'))
    require(local_test and Path(value).is_absolute(),
            'Local file paths require --local-test and must be absolute: ' + value)
    return None, str(Path(value).resolve())


def _source_key(value):
    """Canonicalize redirects and EOS /store aliases for overlap checks."""
    parsed = urlsplit(value)
    if parsed.scheme == 'root':
        path = posixpath.normpath('/' + parsed.path.lstrip('/'))
        if '/store/' in path:
            return 'lfn:' + path[path.index('/store/'):]
        return 'root://' + parsed.netloc.lower() + '/' + path
    return str(Path(value).resolve())


def _command(argv):
    subprocess.run([str(value) for value in argv], check=True)


def _stage(url, destination, local_test=False):
    host, path = _location(url, local_test)
    require(not destination.exists(), 'Staging destination already exists: ' + str(destination))
    if host:
        _command(['xrdcp', '--nopbar', '--retry', '3', '--cksum', 'adler32', url, destination])
    else:
        shutil.copy2(path, destination)
    require(destination.is_file(), 'Copy did not create an input file: ' + url)


def _check_publication_target(url, local_test=False):
    host, path = _location(url, local_test)
    if not host:
        require(not Path(path).exists(), 'Publication destination already exists: ' + url)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        return
    result = subprocess.run(['xrdfs', host, 'stat', path], capture_output=True, text=True)
    require(result.returncode != 0, 'Publication destination already exists: ' + url)
    diagnostic = result.stdout + result.stderr
    require('No such file or directory' in diagnostic or '[3011]' in diagnostic,
            'Could not establish that publication destination is absent: ' + diagnostic.strip())
    _command(['xrdfs', host, 'mkdir', '-p', posixpath.dirname(path)])


def _publish(source, url, local_test=False):
    host, path = _location(url, local_test)
    if host:
        # No --force: an existing result is never replaced. POSC removes an
        # incomplete remote file when the transfer itself fails.
        _command(['xrdcp', '--nopbar', '--retry', '3', '--cksum', 'adler32',
                  '--posc', source, url])
    else:
        with Path(source).open('rb') as incoming, Path(path).open('xb') as outgoing:
            shutil.copyfileobj(incoming, outgoing)


def _raw_integrity(path):
    if path.stat().st_size == 0:
        return 'zero-byte ROOT file'
    ROOT = root_merge._root()
    file = None
    try:
        try:
            file = ROOT.TFile.Open(str(path), 'READ')
        except OSError as error:
            return 'unreadable ROOT file: ' + str(error)
        if not file or not file.IsOpen():
            return 'unreadable ROOT file'
        if file.IsZombie():
            return 'zombie ROOT file'
        if file.TestBit(ROOT.TFile.kRecovered):
            return 'recovered ROOT file'
        return None
    finally:
        if file:
            file.Close()


def _check_data_card(path, data_type):
    # Integrity was checked separately: missing or incompatible cards must fail,
    # and may never enter the corrupt-input exclusion policy.
    ROOT = root_merge._root()
    file = ROOT.TFile.Open(str(path), 'READ')
    try:
        require(bool(file) and not file.IsZombie() and not file.TestBit(ROOT.TFile.kRecovered),
                'ROOT input changed after integrity validation: ' + str(path))
        card = file.GetDirectory('JCard')
        require(bool(card), 'Missing JCard: ' + str(path))
        for key, expected in (('DataType', data_type), ('BackgroundMethods', 6)):
            value = card.Get(key)
            require(bool(value), 'Missing JCard/' + key + ': ' + str(path))
            require(root_merge._card_value(value).get('value') == [expected],
                    'Unexpected JCard/' + key + ': ' + str(path))
    finally:
        if file:
            file.Close()



def _promote_event_counter(path):
    """Preserve exact event counts when a data reduction exceeds float precision.

    This modifies only the staged scratch copy, after parent hash verification.
    Every input is normalized before merge_root compares histogram schemas.
    """
    ROOT = root_merge._root()
    file = ROOT.TFile.Open(str(path), 'UPDATE')
    try:
        require(bool(file) and not file.IsZombie() and not file.TestBit(ROOT.TFile.kRecovered),
                'Cannot normalize staged event counter: ' + str(path))
        events = file.Get('nEvents')
        require(bool(events) and events.InheritsFrom('TH1') and events.GetDimension() == 1,
                'Missing or incompatible nEvents histogram: ' + str(path))
        if events.ClassName() == 'TH1D':
            return False
        require(events.ClassName() == 'TH1F', 'Unsupported nEvents storage class: ' + str(events.ClassName()))
        bins = int(events.GetNbinsX())
        promoted = ROOT.TH1D('__eec_double_event_counter', events.GetTitle(), bins,
                             events.GetXaxis().GetXmin(), events.GetXaxis().GetXmax())
        promoted.SetDirectory(0)
        for source_axis, target_axis in ((events.GetXaxis(), promoted.GetXaxis()),
                                         (events.GetYaxis(), promoted.GetYaxis()),
                                         (events.GetZaxis(), promoted.GetZaxis())):
            source_axis.Copy(target_axis)
        if events.GetSumw2N():
            promoted.Sumw2()
        for index in range(events.GetNcells()):
            promoted.SetBinContent(index, events.GetBinContent(index))
            if events.GetSumw2N():
                promoted.SetBinError(index, events.GetBinError(index))
        statistics = array('d', [0.] * 4)
        events.GetStats(statistics)
        promoted.PutStats(statistics)
        promoted.SetEntries(events.GetEntries())
        promoted.SetName('nEvents')
        file.cd()
        require(promoted.Write('nEvents', ROOT.TObject.kOverwrite) > 0,
                'Could not store double-precision event counter')
        return True
    finally:
        if file:
            file.Close()


def _validate_manifest(manifest, local_test):
    require(isinstance(manifest, dict), 'Manifest must be an object')
    require(isinstance(manifest.get('node_id'), str) and
            re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', manifest['node_id']), 'Invalid node_id')
    require(manifest.get('sample') in SAMPLE_TYPES, 'Unsupported data sample')
    require(manifest.get('data_type') == SAMPLE_TYPES[manifest['sample']], 'Sample/DataType mismatch')
    require(manifest.get('kind') in ('batch', 'part', 'sample'), 'Unsupported merge node kind')
    inputs = manifest.get('inputs')
    require(isinstance(inputs, list) and bool(inputs), 'Empty or invalid input list')
    seen = set()
    parent_flags = []
    for entry in inputs:
        require(isinstance(entry, dict), 'Invalid input entry')
        _location(entry.get('url'), local_test)
        key = _source_key(entry['url'])
        require(key not in seen, 'Duplicate immediate input: ' + entry['url'])
        seen.add(key)
        size = entry.get('size_bytes')
        is_parent = 'report_url' in entry
        parent_flags.append(is_parent)
        require((is_parent and size is None) or
                (type(size) is int and size >= 0), 'Invalid input size: ' + entry['url'])
        if is_parent:
            _location(entry.get('report_url'), local_test)
    require(all(parent_flags) or not any(parent_flags), 'Raw and parent inputs cannot be mixed in one node')
    for key in ('output_url', 'report_url'):
        _location(manifest.get(key), local_test)
        require(_source_key(manifest[key]) not in seen, 'Output aliases an input')
    require(_source_key(manifest['output_url']) != _source_key(manifest['report_url']),
            'ROOT and report destinations alias each other')


def _parent_report(path, entry, manifest):
    parent = json.loads(path.read_text())
    require(parent.get('schema_version') == 1 and parent.get('status') in SUCCESS_STATUSES,
            'Parent merge has no successful report: ' + entry['url'])
    require(parent.get('sample') == manifest['sample'] and
            parent.get('data_type') == manifest['data_type'] and
            parent.get('background_methods') == 6, 'Parent data configuration differs')
    publication = parent.get('publication', {})
    require(publication.get('root_published') is True and publication.get('report_published') is True,
            'Parent report does not certify completed publication')
    require(_source_key(publication.get('root_url', '')) == _source_key(entry['url']) and
            _source_key(publication.get('report_url', '')) == _source_key(entry['report_url']),
            'Parent publication URLs differ from the frozen node inputs')
    artifact = parent.get('artifact', {})
    require(type(artifact.get('size_bytes')) is int and artifact['size_bytes'] > 0 and
            isinstance(artifact.get('sha256'), str) and
            re.fullmatch(r'[0-9a-f]{64}', artifact['sha256']), 'Invalid parent artifact metadata')
    included = parent.get('included_source_urls')
    excluded = parent.get('excluded_files')
    require(isinstance(included, list) and bool(included) and
            parent.get('included_source_count') == len(included), 'Invalid parent original-source membership')
    require(isinstance(excluded, list), 'Invalid parent exclusion list')
    require((parent['status'] == 'PARTIAL') == bool(excluded), 'Parent status/exclusion mismatch')
    for value in included:
        require(isinstance(value, str) and bool(value), 'Invalid parent included URL')
    for row in excluded:
        require(isinstance(row, dict) and isinstance(row.get('url'), str) and bool(row['url']) and
                isinstance(row.get('reason'), str) and bool(row['reason']), 'Invalid parent exclusion record')
    keys = [_source_key(value) for value in included] + [_source_key(row['url']) for row in excluded]
    require(len(keys) == len(set(keys)), 'Parent repeats or both includes and excludes an original source')
    totals = parent.get('root_totals', {})
    for key in TOTAL_KEYS:
        value = totals.get(key)
        require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
                'Invalid parent ROOT total: ' + key)
    bins = totals.get('n_events_bins')
    require(isinstance(bins, list) and len(bins) >= 11 and
            all(type(value) in (int, float) and math.isfinite(value) for value in bins),
            'Invalid parent nEvents bins')
    return parent


def _totals_from_parents(parents):
    if not parents:
        return None
    totals = {key: math.fsum(parent['root_totals'][key] for parent in parents) for key in TOTAL_KEYS}
    counts = {len(parent['root_totals']['n_events_bins']) for parent in parents}
    require(len(counts) == 1, 'Parent nEvents bin counts differ')
    totals['n_events_bins'] = [math.fsum(values) for values in
                              zip(*(parent['root_totals']['n_events_bins'] for parent in parents))]
    return totals


def _write_report(work, report):
    (work / 'good_files.txt').write_text(''.join(url + '\n' for url in report['included_source_urls']))
    (work / 'bad_files.txt').write_text(''.join(row['url'] + '\n' for row in report['excluded_files']))
    (work / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


def run(manifest_path, work_dir, *, local_test=False, publish=True):
    work = Path(work_dir).resolve()
    require(local_test or not str(work).startswith(('/afs/', '/eos/')), 'Worker requires local scratch')
    work.mkdir(parents=True, exist_ok=True)
    require(not any(work.iterdir()), 'Worker scratch directory must be empty: ' + str(work))
    output = work / 'result.root'
    report = {'schema_version': 1, 'node_id': None, 'sample': None, 'kind': None,
              'data_type': None, 'background_methods': 6, 'status': 'FAIL',
              'snapshot': {}, 'included_source_urls': [], 'included_source_count': 0,
              'excluded_files': [], 'immediate_sources': [], 'artifact': None,
              'root_totals': None, 'event_counter_normalization': {
                  'histogram': 'nEvents', 'target_class': 'TH1D', 'promoted_input_count': 0,
                  'scope': 'Staged copies only; preserve bins, stored errors, axes, entries and statistics before exact-count merging'},
              'publication': {'root_url': None, 'report_url': None,
                                                 'root_published': False, 'report_published': False}}
    try:
        manifest_path = Path(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        for key in ('node_id', 'sample', 'kind', 'data_type', 'snapshot'):
            if key in manifest:
                report[key] = manifest[key]
        report['manifest_sha256'] = sha256(manifest_path)
        report['publication'].update(root_url=manifest.get('output_url'), report_url=manifest.get('report_url'))
        _validate_manifest(manifest, local_test)
        local_inputs, parents, original_keys = [], [], set()
        for index, entry in enumerate(manifest['inputs']):
            immediate = dict(entry)
            immediate['status'] = 'PENDING'
            report['immediate_sources'].append(immediate)
            parent = None
            if 'report_url' in entry:
                parent_path = work / ('parent_%06d.json' % index)
                _stage(entry['report_url'], parent_path, local_test)
                parent = _parent_report(parent_path, entry, manifest)
                original = parent['included_source_urls'] + [row['url'] for row in parent['excluded_files']]
                keys = {_source_key(value) for value in original}
                require(not original_keys.intersection(keys), 'Parent parts overlap in original source URLs')
                original_keys.update(keys)
                immediate['parent_node_id'] = parent.get('node_id')
                immediate['parent_report_sha256'] = sha256(parent_path)
            local = work / ('input_%06d.root' % index)
            if parent:
                _stage(entry['url'], local, local_test)
            else:
                immediate['source_metadata'] = stage_frozen_input(
                    entry['url'], local, entry,
                    lambda: _stage(entry['url'], local, local_test), local_test)
            actual_size = local.stat().st_size
            if entry.get('size_bytes') is not None:
                require(actual_size == entry['size_bytes'], 'Input size changed since snapshot: ' + entry['url'])
            immediate['staged_size_bytes'] = actual_size
            if parent:
                require(actual_size == parent['artifact']['size_bytes'] and
                        sha256(local) == parent['artifact']['sha256'], 'Parent artifact size/SHA256 mismatch')
                reason = _raw_integrity(local)
                require(reason is None, 'Previously merged parent ROOT is invalid: ' + str(reason))
            else:
                reason = _raw_integrity(local)
                if reason is not None:
                    immediate.update(status='EXCLUDED', reason=reason)
                    report['excluded_files'].append({'url': entry['url'], 'reason': reason,
                                                     **{key: entry[key] for key in ('job_id', 'request_name') if key in entry}})
                    print('MERGE_RAW_EXCLUDED ' + json.dumps(report['excluded_files'][-1]), flush=True)
                    continue
            _check_data_card(local, manifest['data_type'])
            promoted = _promote_event_counter(local)
            immediate['event_counter_promoted'] = promoted
            report['event_counter_normalization']['promoted_input_count'] += int(promoted)
            local_inputs.append(local)
            immediate['status'] = 'INCLUDED'
            if parent:
                parents.append(parent)
                report['included_source_urls'].extend(parent['included_source_urls'])
                report['excluded_files'].extend(parent['excluded_files'])
            else:
                report['included_source_urls'].append(entry['url'])
            report['included_source_count'] = len(report['included_source_urls'])
            print('MERGE_INPUT_VERIFIED %d/%d' % (index + 1, len(manifest['inputs'])), flush=True)
        require(bool(local_inputs), 'No valid ROOT inputs remain after raw integrity exclusions')
        validation = root_merge.merge_root(local_inputs, output, _totals_from_parents(parents))
        report['root_totals'] = validation['expected_totals']
        report['validation'] = {key: validation[key] for key in
                                ('status', 'input_count', 'validation_scope', 'float_bin_relative_tolerance')}
        report['artifact'] = {'sha256': sha256(output), 'size_bytes': output.stat().st_size}
        report['status'] = 'PARTIAL' if report['excluded_files'] else 'PASS'
        report['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        if publish:
            # Check both destinations before publishing the ROOT artifact.
            _check_publication_target(manifest['output_url'], local_test)
            _check_publication_target(manifest['report_url'], local_test)
            _publish(output, manifest['output_url'], local_test)
            report['publication']['root_published'] = True
            report['publication']['report_published'] = True
            _write_report(work, report)
            try:
                _publish(work / 'result.json', manifest['report_url'], local_test)
            except Exception:
                report['publication']['report_published'] = False
                raise
        print('MERGE_NODE_SUCCESS ' + json.dumps({'node_id': report['node_id'], 'status': report['status'],
                                                'included': report['included_source_count'],
                                                'excluded': len(report['excluded_files'])}), flush=True)
        return report
    except Exception as error:
        report['status'] = 'FAIL'
        report['error'] = {'type': type(error).__name__, 'message': str(error)}
        raise
    finally:
        report['included_source_count'] = len(report['included_source_urls'])
        report['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        # A zero-byte placeholder on failure is explicitly not a valid artifact;
        # it permits Condor to return the same diagnostic filenames every time.
        if not output.exists():
            output.touch()
        _write_report(work, report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--local-test', action='store_true', help='Allow absolute local paths for test fixtures')
    parser.add_argument('--no-publish', action='store_true')
    args = parser.parse_args()
    try:
        run(args.manifest, args.work_dir, local_test=args.local_test, publish=not args.no_publish)
    except Exception as error:
        print('MERGE_NODE_FAILED ' + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
