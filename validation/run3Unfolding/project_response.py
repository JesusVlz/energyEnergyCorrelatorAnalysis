#!/usr/bin/env python3
"""Project existing response histograms with the standard macro and fresh libraries.

This performs projection only. It does not run the analyzer or unfolding.
A new scratch directory avoids replacing plotting/DrawingClasses.so, whose
external dependencies may differ from the installed ROOT environment.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import tempfile
import subprocess

REPO = Path(__file__).resolve().parents[2]


def project(input_path, output_path, power, mask, scratch):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    input_path = input_path.resolve(strict=True)
    output_path = output_path.resolve()
    if output_path.exists():
        raise ValueError('Use a new projection output; the standard macro writes UPDATE: '+str(output_path))
    scratch = scratch.resolve()
    scratch.mkdir(parents=True, exist_ok=False)
    sources = [REPO/'plotting'/(name+'.cxx') for name in (
        'EECHistogramManager', 'EECCard', 'AlgorithmLibrary', 'EECBackgroundScale',
        'EECSignalToBackgroundUnfoldingScale', 'EECUnfoldConfiguration')]
    sources += [REPO/'src/ConfigurationCard.cxx', REPO/'src/EECHistograms.cxx']
    library = scratch/'standardManager.so'
    command = ['g++', '-shared', '-fPIC', '-O0', '-I'+str(REPO),
               *map(str, sources),
               *shlex.split(subprocess.check_output(['root-config', '--cflags', '--libs'], text=True)),
               '-lEG', '-o', str(library)]
    with (scratch/'build.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    if ROOT.gSystem.Load(str(library)) < 0:
        raise RuntimeError('Cannot load freshly compiled standard manager')
    macro_path = REPO/'plotting/projectEEChistograms.C'
    original = macro_path.read_text()
    directive = 'R__LOAD_LIBRARY(plotting/DrawingClasses.so)'
    if original.count(directive) != 1:
        raise ValueError('Unexpected standard projection macro library directive')
    ROOT.gInterpreter.AddIncludePath(str(REPO/'plotting'))
    if not ROOT.gInterpreter.Declare(original.replace(directive, '')):
        raise RuntimeError('Cannot declare standard projection macro')
    source = ROOT.TFile.Open(str(input_path), 'READ')
    if not source or source.IsZombie():
        raise ValueError('Cannot open projection input')
    try:
        required = ['JCard']
        if mask & 8192:
            required += ['jetPtUnfoldingMeasured', 'jetPtUnfoldingTruth', 'jetPtUnfoldingResponse']
        if mask & 32768:
            required += ['oneDimensionalJetPtUnfoldingMeasured', 'oneDimensionalJetPtUnfoldingTruth',
                         'oneDimensionalJetPtUnfoldingResponse']
        if mask & 65536:
            required += ['jetPtUnfoldingCovariance']
        for name in required:
            if not source.Get(name):
                raise ValueError('Missing required raw projection object: '+name)
    finally:
        source.Close()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.TemporaryDirectory(prefix='eec-response-projection-')
    staged_output = Path(temporary.name)/output_path.name
    ROOT.projectEEChistograms(str(input_path), str(staged_output), mask, power)
    projected = ROOT.TFile.Open(str(staged_output), 'READ')
    if not projected or projected.IsZombie():
        raise RuntimeError('Standard macro did not write a valid projection output')
    try:
        for name in required:
            if not projected.Get(name):
                raise RuntimeError('Projected output is missing '+name)
    finally:
        projected.Close()
    receipt = dict(input=str(input_path), output=str(output_path), power=power, mask=mask,
                   root_version=ROOT.gROOT.GetVersion(), build_command=command,
                   macro=str(macro_path), macro_sha256=hashlib.sha256(original.encode()).hexdigest(),
                   source_sha256={str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [*sources, Path(__file__).resolve(), macro_path]},
                   macro_change='Removed only stale DrawingClasses library load directive; function body unchanged.',
                   scope='Projection only; no analyzer, response training, unfolding or iteration choice.')
    # EECCard's known-key list omits selection metadata. Preserve extra keys
    # without replacing the projected known keys (notably WeightExponent).
    source = ROOT.TFile.Open(str(input_path), 'READ')
    projected = ROOT.TFile.Open(str(staged_output), 'UPDATE')
    added, original_names = [], []
    try:
        original_card = source.GetDirectory('JCard')
        output_card = projected.GetDirectory('JCard')
        input_card = projected.mkdir('InputJCard')
        if not original_card or not output_card or not input_card:
            raise RuntimeError('Cannot preserve complete input JCard metadata')
        for key in original_card.GetListOfKeys():
            name = key.GetName()
            original_names.append(name)
            obj = key.ReadObj()
            if not output_card.GetKey(name):
                output_card.cd()
                obj.Write(name)
                added.append(name)
            input_card.cd()
            obj.Write(name)
        receipt['extra_jcard_keys_preserved'] = added
        receipt['original_card_copy'] = 'InputJCard'
        checksum = hashlib.sha256()
        with input_path.open('rb') as stream:
            for block in iter(lambda: stream.read(8*1024*1024), b''):
                checksum.update(block)
        receipt['input_sha256'] = checksum.hexdigest()
        projected.cd()
        ROOT.TObjString(json.dumps(receipt, sort_keys=True)).Write('ProjectionProvenance')
    finally:
        projected.Close()
        source.Close()
    # Keep repeated ROOT UPDATE operations local. Publish the completed file in
    # one copy, rather than update a remote EOS file through several handles.
    verified = ROOT.TFile.Open(str(staged_output), 'READ')
    if not verified or verified.IsZombie() or verified.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Staged projection requires recovery or cannot be reopened')
    try:
        for name in original_names:
            if not verified.GetDirectory('InputJCard').GetKey(name):
                raise RuntimeError('Staged original card lost '+name)
            if not verified.GetDirectory('JCard').GetKey(name):
                raise RuntimeError('Staged projected card lost '+name)
        if not verified.GetKey('ProjectionProvenance'):
            raise RuntimeError('Staged projection lost its provenance')
    finally:
        verified.Close()
    staged_hash = hashlib.sha256(staged_output.read_bytes()).hexdigest()
    shutil.copyfile(staged_output, output_path)
    if hashlib.sha256(output_path.read_bytes()).hexdigest() != staged_hash:
        raise RuntimeError('Published projection SHA256 differs from validated staged file')
    receipt['output_sha256'] = staged_hash
    receipt['publication'] = 'Validated local ROOT staging, then one destination copy with SHA256 verification'
    temporary.cleanup()
    (scratch/'projection_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'PROJECTED', 'output': str(output_path), 'power': power, 'mask': mask}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--power', type=int, choices=(1, 2), required=True)
    parser.add_argument('--mask', type=int, default=108803)
    parser.add_argument('--scratch', type=Path, required=True, help='new directory for isolated build and receipt')
    args = parser.parse_args()
    if args.mask <= 0 or args.mask >= 1 << 17:
        parser.error('--mask must select existing standard projection bits0..16')
    project(args.input, args.output, args.power, args.mask, args.scratch)


if __name__ == '__main__':
    main()
