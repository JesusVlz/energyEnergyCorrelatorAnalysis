#!/usr/bin/env python3
"""One Condor pair task; native exports exist only in worker scratch.

The parent campaign submits and stages files. This worker never submits jobs.
Output consists of compact counts/cache/identity shards and their receipt.
"""
import argparse,json,shutil,sys,tempfile,time
from collections import deque
from pathlib import Path
from types import SimpleNamespace
import batch_pairs as batch
import export_forest
import pair_counts as p


def report_native_failure(work,prefix,error):
    """Retain native diagnostics before TemporaryDirectory removes the work area."""
    print('PAIR_NATIVE_FAILURE: '+repr(error),file=sys.stderr,flush=True)
    for log in sorted(Path(work).glob('*.native/*.log')):
        saved=Path(str(prefix)+'.'+log.stem+'.native_failure.log')
        try:
            shutil.copyfile(log,saved)
            with log.open(errors='replace') as stream:tail=list(deque(stream,maxlen=80))
            print('Native log: '+str(log)+'; retained copy: '+str(saved),file=sys.stderr)
            print(''.join(tail),file=sys.stderr,end='' if tail and tail[-1].endswith('\n') else '\n',flush=True)
        except OSError as log_error:
            print('Could not retain native log '+str(log)+': '+str(log_error),file=sys.stderr,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--task-index',type=int,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--scratch',type=Path,required=True)
    parser.add_argument('--input-local',type=Path,help='already staged byte-identical input; original identity remains in task provenance')
    args=parser.parse_args();started=time.monotonic();manifest=batch.read_json(args.manifest)
    if not 0<=args.task_index<len(manifest['tasks']):parser.error('task index outside frozen manifest')
    task=manifest['tasks'][args.task_index];original=Path(task['input']);source=args.input_local or original
    identity=task.get('input_identity',{})
    if source.stat().st_size!=identity.get('size_bytes',source.stat().st_size):raise ValueError('Frozen input size mismatch')
    if args.input_local is None and identity.get('mtime_ns',source.stat().st_mtime_ns)!=source.stat().st_mtime_ns:raise ValueError('Frozen input modification time mismatch')
    args.scratch.mkdir(parents=True,exist_ok=True);args.output_dir.mkdir(parents=True,exist_ok=True)
    prefix=args.output_dir/f"{task['system']}_{task['file_index']:03d}"
    if Path(str(prefix)+'.receipt.json').exists():raise ValueError('Output task is already complete; campaign must verify before resume')
    with tempfile.TemporaryDirectory(prefix=prefix.name+'-',dir=args.scratch) as temporary:
        output=Path(temporary)/'events.jsonl'
        try:
            export_forest.export(SimpleNamespace(input=source,output=output,system=task['system'],sample_namespace=task['sample_namespace'],
                                                 payload=None,card=None,max_events=task['entry_limit'],responses=False,parity=False,reuse_native=False))
        except Exception as error:
            report_native_failure(temporary,prefix,error)
            raise
        metadata_path=output.with_suffix('.jsonl.metadata.json');metadata=batch.read_json(metadata_path)
        metadata.update(original_input=task['input'],original_input_identity=identity,raw_input_sha256=p.file_digest(source),
                        native_binary_sha256=p.file_digest(export_forest.REPO/'eecAnalysis'),batch_worker_sha256=p.file_digest(Path(__file__)),
                        batch_reducer_sha256=p.file_digest(Path(batch.__file__)),frozen_task_manifest_sha256=p.file_digest(args.manifest))
        batch.write_json(metadata_path,metadata)
        ledger=Path(metadata['native_outputs'])/'reco.jsonl.events.jsonl'
        receipt=batch.train(output,prefix,ledger,metadata_path,task)
    print(json.dumps({'status':'PAIR_TRAINING_SHARD_COMPLETE','prefix':str(prefix),'event_counts':receipt['event_counts'],
                      'persisted_bytes':sum(v['bytes'] for v in receipt['files'].values()),'elapsed_seconds':time.monotonic()-started}))
if __name__=='__main__':main()
