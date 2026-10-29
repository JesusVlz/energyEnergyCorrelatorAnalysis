#!/usr/bin/env python3
"""Build and run synthetic Run-3 ID/veto regressions; no real map is selected."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import time

REPO = Path(__file__).resolve().parents[2]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    destination = args.output.resolve()
    records = []
    def run(command, label):
        started = time.time()
        result = subprocess.run(command, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = destination / (label + '.log')
        log.write_text(result.stdout)
        records.append({'command': command, 'returncode': result.returncode, 'seconds': time.time()-started, 'log': str(log)})
        (destination/'commands.json').write_text(json.dumps(records,indent=2)+'\n')
        print(label + ': ' + ('PASS' if result.returncode == 0 else 'FAIL'), flush=True)
        if result.returncode:
            raise RuntimeError(result.stdout[-6000:])
    run(['make','-j4'], 'build')
    flags = shlex.split(subprocess.check_output(['root-config','--cflags','--libs'],text=True))
    objects = [str(p.relative_to(REPO)) for p in sorted((REPO/'src').glob('*.o'))]
    fixtures = [
        ('core','validation/run3Jets/testRun3JetSelection.cxx',['src/Run3JetSelection.o'],False),
        ('readers','validation/run3Jets/testRun3JetSelectionReaders.cxx',['src/ForestReader.o','src/HighForestReader.o'],False),
        ('analyzer','validation/run3Jets/checkAnalyzerJetSelection.cxx',objects,True),
        ('calibration','validation/run3Jets/checkCalibrationInputs.cxx',objects,True),
        ('response_regression','validation/run3Response/checkAnalyzerConfiguration.cxx',objects,True),
    ]
    for label,source,inputs,private in fixtures:
        binary = destination / label
        command = ['g++','-g','-Wall','-Isrc','-I.']
        if private:
            command += ['-fno-access-control']
        command += [source,*inputs,*flags,'-lEG','-o',str(binary)]
        run(command,label+'_compile')
        run([str(binary),str(destination/(label+'_fixtures'))],label)
    run(['python3','-m','unittest','discover','-s','validation/run3Jets','-p','test_*.py','-v'],'python_diagnostics')
    print('PASS: synthetic core, reader, analyzer, response and diagnostic checks')

if __name__ == '__main__':
    main()
