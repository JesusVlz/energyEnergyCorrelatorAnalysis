#!/usr/bin/env python3
"""Install the three original ppRef ROOT tables; no network, decoding, or job submission.

Default: dry run. --write copies verified bytes. --check verifies installed files.
--activate-base-cards also enables nominal pair correction in the six pp R4
compatibility cards, after validating their configuration and all three payloads.
It never changes OO cards, response guards, calibration files, or ROOT contents.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tarfile

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / 'trackCorrectionTables/pp2024/pair/20261116/manifest.json'
BASE_CARDS = (
    'cardEECpp5v36TeV.input', 'cardEECpp5v36TeVMC.input',
    'cardEECpp5v36TeVR04.input', 'cardEECpp5v36TeVMCR04.input',
    'crab/cardEECpp5v36TeV.input', 'crab/cardEECpp5v36TeVMC.input',
)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def manifest():
    result = json.loads(MANIFEST.read_text())
    if result['schema'] != 'ppref-pair-install-v1':
        raise ValueError('Unsupported installation manifest')
    if set(result['tables']) != {'Nominal', 'Loose', 'Tight'}:
        raise ValueError('Require all three working points')
    for row in result['tables'].values():
        for key in ('path', 'single_track_path'):
            p = Path(row[key])
            if p.is_absolute() or '..' in p.parts:
                raise ValueError('Unsafe manifest path')
    return result


def verify_bytes(data, row):
    if len(data) != row['bytes'] or digest(data) != row['sha256'] or data[:4] != b'root':
        raise ValueError('Payload bytes/checksum do not match: ' + row['path'])


def verify_single_track(row):
    path = REPO / row['single_track_path']
    if not path.is_file() or digest(path.read_bytes()) != row['single_track_sha256']:
        raise ValueError('Single-track payload differs from pair derivation: ' + str(path))


def verify_root(path, wp, row, spec):
    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError('Use the analysis ROOT/PyROOT environment for metadata checks') from exc
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie():
        raise ValueError('Unreadable ROOT table: ' + str(path))
    try:
        tags = {'Run3PairCorrectionVersion': spec['pair_version'],
                'Run3TrackSelectionVersion': spec['track_selection_version'],
                'Run3TrackSelection': wp, 'Run3PairCollisionSystem': 'pp',
                'Run3SingleTrackPayloadSHA256': row['single_track_sha256']}
        for key, expected in tags.items():
            obj = f.Get(key)
            if not obj or str(obj.GetTitle()) != expected:
                raise ValueError('Missing/mismatched ROOT tag ' + key)
        for key, expected in row['axes'].items():
            obj = f.Get('JCard/' + key)
            if not obj or obj.GetNrows() != len(expected):
                raise ValueError('Missing/mismatched ROOT axis ' + key)
            actual = [float(obj.GetMatrixArray()[i]) for i in range(obj.GetNrows())]
            if actual != expected:
                raise ValueError('Different ROOT axis ' + key)
        for key in f.GetListOfKeys():
            name = key.GetName()
            if not name.startswith('trackPairEfficiencyCorrection') or not name.endswith('_valid'):
                continue
            mask, hist = f.Get(name), f.Get(name[:-6])
            if not hist or hist.GetNbinsX() != mask.GetNbinsX():
                raise ValueError('Missing correction histogram for ' + name)
            for i in range(1, hist.GetNbinsX()+1):
                value = hist.GetBinContent(i)
                positive = math.isfinite(value) and value > 0
                if bool(mask.GetBinContent(i)) != positive:
                    raise ValueError('Validity mask/content mismatch: ' + name)
    finally:
        f.Close()


def card_values(text):
    result = {}
    for line in text.splitlines():
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        if len(fields) < 2 or fields[0] in result:
            raise ValueError('Malformed/duplicate card key: ' + fields[0])
        result[fields[0]] = ' '.join(fields[1:])
    return result


def activated_card(text, payload):
    c = card_values(text)
    if int(c['DataType']) not in (11, 12) or float(c['JetRadius']) != .4:
        raise ValueError('Activation is ppRef R4 only')
    if int(c['McCorrelationType']) != 0 or int(c['FilledHistograms']) & (32 | 64 | 128 | 256):
        raise ValueError('Use an ordinary Data/RecoReco EEC card, not responses or systematics')
    if c.get('Run3TrackSelection', 'Nominal') != 'Nominal' or float(c['VertexMaxDistance']) != 3:
        raise ValueError('Compatibility-card activation is Nominal only')
    if float(c['MaxTrackPtRelativeError']) != .1 or float(c['MaxTrackPtCut']) > 500:
        raise ValueError('Track cuts do not match the supplied map')
    if float(c['MinJetPtCut']) < 5 or float(c['MaxJetPtCut']) > 5020:
        raise ValueError('Jet selection exceeds the stored map axes')
    changes = {'Run3TrackPairCorrectionFile': payload, 'DisableTrackPairEfficiencyCorrection': '0',
               'Run3JetID': '1', 'Run3UsePpL1FastJet': '0', 'MinPairTrackPt': '2'}
    lines = []
    for line in text.splitlines():
        fields = line.split('#', 1)[0].split()
        if fields and fields[0] in changes:
            lines.append(fields[0] + ' ' + changes.pop(fields[0]))
        else:
            lines.append(line)
    lines += [k + ' ' + v for k, v in changes.items()]
    warning = '# PRELIMINARY pp pair table 20261116: original derivation JetID OFF; JetID-on closure pending.'
    if warning not in lines:
        lines.insert(0, warning)
    return '\n'.join(lines) + '\n'


def atomic_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.install-tmp')
    with tmp.open('xb') as stream:
        stream.write(data)
    os.replace(tmp, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, help='Original uploaded tracking-closure tar.gz')
    parser.add_argument('--write', action='store_true', help='Install checked ROOT bytes')
    parser.add_argument('--check', action='store_true', help='Verify installed files and ROOT metadata')
    parser.add_argument('--activate-base-cards', action='store_true', help='Also enable the six nominal pp compatibility cards')
    args = parser.parse_args()
    if args.activate_base_cards and not args.write:
        parser.error('--activate-base-cards requires --write')
    if not (args.archive or args.check):
        parser.error('Supply --archive or --check')
    spec = manifest()
    pending = {}
    if args.archive:
        with tarfile.open(args.archive, 'r:gz') as archive:
            for wp, row in spec['tables'].items():
                member = archive.getmember(row['source_member'])
                if not member.isfile() or member.size != row['bytes']:
                    raise ValueError('Wrong archive member: ' + member.name)
                data = archive.extractfile(member).read()
                verify_bytes(data, row)
                verify_single_track(row)
                path = REPO / row['path']
                if path.exists() and path.read_bytes() != data:
                    raise FileExistsError('Refuse to overwrite a different payload: ' + str(path))
                pending[path] = data
    card_changes = {}
    if args.activate_base_cards:
        for name in BASE_CARDS:
            path = REPO / name
            card_changes[path] = activated_card(path.read_text(), spec['tables']['Nominal']['path'])
    if args.write:
        for path, data in pending.items():
            if not path.exists():
                atomic_write(path, data)
    if args.check or args.write:
        for wp, row in spec['tables'].items():
            path = REPO / row['path']
            if not path.is_file():
                raise FileNotFoundError('Install the original ROOT payload first: ' + str(path))
            verify_bytes(path.read_bytes(), row)
            verify_single_track(row)
            verify_root(path, wp, row, spec)
        for path, text in card_changes.items():
            if path.read_text() != text:
                atomic_write(path, text.encode())
    print(json.dumps({'status': 'VERIFIED_LOCAL_INSTALL' if args.write or args.check else 'DRY_RUN',
                      'files': [r['path'] for r in spec['tables'].values()],
                      'activated_cards': [str(p.relative_to(REPO)) for p in card_changes],
                      'jet_id': 1, 'physics_closure': 'PENDING_JETID_ON',
                      'git_commit_or_push': 'NOT_PERFORMED'}, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, RuntimeError, tarfile.TarError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        sys.exit(1)
