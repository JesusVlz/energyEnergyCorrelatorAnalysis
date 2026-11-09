#!/usr/bin/env python3
"""Verify the Run-3 JEC payloads used by the EEC analysis."""
import hashlib
from pathlib import Path
import sys

PAYLOADS = {
    'OO2025/Prompt25HIOO_V1_MC_L2Relative_AK4PF.txt':
        'd3f8498ad7ee1bc5de1dcf0de802b5d22d0ce5d5',
    'pp2024/Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt':
        '7bcc3847f0f5e7d3b0b927e6df64d2d7c9ff60fa',
    'pp2024/Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt':
        '2f2f4b04cb1134919b3a54e6ba5d11a73e3ac710',
    'pp2024/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt':
        '99ea207b9977fd96fe40b4833f81d3fa5883ce8a',
}

def git_blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()

def main():
    root = Path(__file__).resolve().parents[2] / 'jetEnergyCorrections'
    for relative, expected in PAYLOADS.items():
        data = (root / relative).read_bytes()
        actual = git_blob(data)
        if actual != expected:
            raise ValueError(f'{relative}: git blob {actual}, expected {expected}')
        first = data.decode('utf-8').splitlines()[0]
        if 'JetEta' not in first or 'JetPt' not in first:
            raise ValueError(f'{relative}: expected JetEta/JetPt JEC payload')
        if 'L1FastJet' in relative and ('Rho' not in first or 'JetA' not in first):
            raise ValueError(f'{relative}: expected Rho and JetA dependencies')
        print(f'PASS {relative}: git_blob={actual} sha256={hashlib.sha256(data).hexdigest()}')
    print('PASS: ppRef payload set contains optional L1FastJet, active L2Relative V2, and data residual.')

if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
