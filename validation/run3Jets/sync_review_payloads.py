#!/usr/bin/env python3
"""Verify the requested Run-3 payloads; --fetch retrieves pinned upstream bytes.

This does NOT change the analyzer's active correction chain. Run from any directory.
All downloads are verified before any local file is changed. No ROOT dependency.
"""
import argparse
import hashlib
from pathlib import Path
import sys
import urllib.request

SOURCE = 'thenicholasbarnett/HIN-Jets'
REVISION = 'c64293a132af388af1827044779af36e57f884d2'
PAYLOADS = {
    'Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt': ('2025OO', 'd3f8498ad7ee1bc5de1dcf0de802b5d22d0ce5d5'),
    'Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt': ('2024ppRef', '7bcc3847f0f5e7d3b0b927e6df64d2d7c9ff60fa'),
    'Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt': ('2024ppRef', '2f2f4b04cb1134919b3a54e6ba5d11a73e3ac710'),
    'Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt': ('2024ppRef', '99ea207b9977fd96fe40b4833f81d3fa5883ce8a'),
}


def git_blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode('ascii') + b'\0' + data).hexdigest()


def verify(name, data, expected):
    actual = git_blob(data)
    if actual != expected:
        raise ValueError(f'{name}: blob {actual}, expected {expected}; refuse to certify this file')
    text = data.decode('utf-8')
    if not text.startswith('{1 JetEta '):
        raise ValueError(f'{name}: unexpected header')
    for number, line in enumerate(text.splitlines()[1:], 2):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 3 or len(fields) != 3 + int(fields[2]):
            raise ValueError(f'{name}:{number}: invalid declared record size')
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch', action='store_true', help='retrieve exact pinned source bytes')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    destination = root / 'jetEnergyCorrections/run3/HIN-Jets_c64293a'
    verified = {}
    for name, (era, expected) in PAYLOADS.items():
        if args.fetch:
            url = f'https://raw.githubusercontent.com/{SOURCE}/{REVISION}/txt/{era}/{name}'
            request = urllib.request.Request(url, headers={'User-Agent': 'OO-EEC-pinned-JEC-import'})
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
        else:
            data = (destination / name).read_bytes()
        digest = verify(name, data, expected)
        verified[name] = data
        print(f'PASS {name}: git_blob={expected} sha256={digest} bytes={len(data)}')
    if args.fetch:
        destination.mkdir(parents=True, exist_ok=True)
        for name, data in verified.items():
            target = destination / name
            temporary = target.with_suffix(target.suffix + '.tmp')
            temporary.write_bytes(data)
            temporary.replace(target)
    print('All four payloads verified. Active analyzer configuration was not modified.')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
