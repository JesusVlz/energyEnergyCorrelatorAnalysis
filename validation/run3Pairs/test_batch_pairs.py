#!/usr/bin/env python3
"""Batch equivalence and global duplicate gates; no physics closure claim."""
import copy,json,tempfile,unittest
from pathlib import Path
from collections import defaultdict
import batch_pairs as b
import pair_counts as p
import plot_closure as plot
from test_pair_counts import event,jet,track

class BatchTest(unittest.TestCase):
    def fixtures(self,directory):
        tasks=[];all_events=[]
        for file in range(3):
            values=[]
            for i in range(10):
                number=file*100+i
                js=[jet(.012,reco=i%3!=0),jet(.2),jet(.06,pt=100.)]
                e=event(number,js);e['event_weight']=1.+i/10;e['source_entry']=i
                e.update(gen=[track(0,3.1,0),track(1,4.1,.012),track(2,3.1,1),track(3,4.1,1.2)],reco=[track(0,3.1,0),track(1,4.1,.012),track(2,3.1,1),track(3,4.1,1.2)])
                e['jets'][0]['reco'] and e['jets'][0]['reco'][0].update(correction=1.25)
                if file==2 and i==9:
                    e['reco'][3].update(correction=0.,unsupported_single=True)
                values.append(e)
            source=directory/f'input{file}.jsonl';source.write_text(''.join(json.dumps(e)+'\n' for e in values))
            ledger=directory/f'ledger{file}.jsonl';ledger.write_text(''.join(json.dumps({'source_entry':e['source_entry'],'id':e['id']})+'\n' for e in values)+json.dumps({'source_entry':10,'id':[2,3,file]})+'\n')
            task={'system':'OO','file_index':file+1,'sample_namespace':'OO-synthetic','input':str(source),'entry_limit':11};tasks.append(task)
            b.train(source,directory/f'OO_{file+1:03d}',ledger,task=task);all_events.extend(values)
        manifest=directory/'manifest.json';manifest.write_text(json.dumps({'tasks':tasks}))
        return manifest,all_events

    def assert_numbers(self,first,second,path=''):
        if isinstance(first,dict):
            self.assertEqual(set(first),set(second),path)
            for k in first:self.assert_numbers(first[k],second[k],path+'/'+str(k))
        elif isinstance(first,list):
            self.assertEqual(len(first),len(second),path)
            for i,(a,c) in enumerate(zip(first,second)):self.assert_numbers(a,c,path+'/'+str(i))
        elif isinstance(first,(float,int)) and not isinstance(first,bool):self.assertAlmostEqual(first,second,delta=1e-10*max(1.,abs(first)),msg=path)
        else:self.assertEqual(first,second,path)

    def test_equivalence_counts_masks_closure_and_covariance(self):
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary);manifest,values=self.fixtures(d);truth=p.run(values,min_pairs=2)
            trained=b.merge_training(manifest,d,'OO',d/'training.json',d/'ids.sqlite',min_pairs=2)
            for name in ['map','jet_independent_map']:
                self.assert_numbers(b.unrows(truth[name]),b.unrows(trained[name]),name)
            for s in trained['shards']:b.validate_cache(s['holdout'],d/'training.json',d/(s['prefix']+'.validation.json.gz'),chunk_size=3)
            final=b.merge_validation(d/'training.json',d,d/'final.json')
            self.assert_numbers([r for r in truth['holdout_spectra'] if r['key'][0]=='OO'],final['holdout_spectra'],'spectra')
            self.assert_numbers(b.unrows(truth['jet_event_sums']),b.unrows(final['jet_event_sums']),'jets')
            self.assert_numbers(b.unrows(truth['correction_application']),b.unrows(final['correction_application']),'use')
            self.assert_numbers(dict(plot.event_statistics(truth)),b.unrows(final['event_statistics']),'absolute covariance')
            _,old_shapes,_=plot.normalized_shapes(truth);new=b.unrows(final['shape_statistics'])
            for key,row in new.items():
                for name,value in row.items():self.assertAlmostEqual(value,old_shapes[key].get(name,0.),delta=1e-10*max(1.,abs(value)),msg=str((key,name)))
            self.assertEqual(final['event_counts']['raw'],33)
            self.assertEqual(final['event_counts']['selected'],30)
            self.assertEqual(final['event_counts']['holdout'],sum(p.fold(v)==1 for v in values))

    def test_duplicate_rejected_raw_event_blocks_global_merge(self):
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary);manifest,_=self.fixtures(d)
            import gzip,sqlite3
            ledger=d/'OO_002.identities.bin.gz'
            with gzip.open(ledger,'rb') as f:records=f.read()
            with gzip.open(ledger,'wb') as f:f.write(records[:-b.IDENTITY.size]+b.IDENTITY.pack(2,3,0,0))
            r=b.read_json(d/'OO_002.receipt.json');r['files']['identities'].update(sha256=p.file_digest(ledger),bytes=ledger.stat().st_size);b.write_json(d/'OO_002.receipt.json',r)
            with self.assertRaises(sqlite3.IntegrityError):b.merge_training(manifest,d,'OO',d/'table.json',d/'ids.sqlite')
            self.assertFalse((d/'table.json').exists())

    def test_changed_cache_blocks_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary);manifest,_=self.fixtures(d);r=b.merge_training(manifest,d,'OO',d/'training.json',d/'ids.sqlite')
            cache=Path(r['shards'][0]['holdout']);cache.write_bytes(cache.read_bytes()+b'changed')
            with self.assertRaisesRegex(ValueError,'bound'):b.validate_cache(cache,d/'training.json',d/'closure.json')

    def test_no_per_file_normalized_corrections_persisted(self):
        with tempfile.TemporaryDirectory() as temporary:
            d=Path(temporary);self.fixtures(d);r=b.read_json(d/'OO_001.counts.json.gz')
            for family in ['differential','jet_independent']:
                for name in ['cells','controls']:
                    for row in r[family][name]:self.assertEqual(set(row),{'key',*b.COUNT_FIELDS})

    def test_sparse_single_event_truth_covariance_uses_checked_variance(self):
        # One physical event fills two angular bins. Its normalized shape has
        # zero statistical variance, but expanded moments cancel below zero.
        channel=('OO','R04',0,1,'MB');spectra=[];event_rows=[]
        for ibin,value in [(3,.03),(4,.3)]:
            key=(*channel,ibin)
            spectra.append({'key':list(key),'status':'MEASURED_HOLDOUT',**{f:value for f in plot.FIELDS}})
            event_rows.append({'key':list(key),**{plot.EVENT_FIELDS.get(f,f):value for f in plot.FIELDS}})
        vectors=[{'jets':[],'channels':event_rows}]
        result={'binning':{'dr_support':list(p.DR)},'holdout_spectra':spectra,'holdout_event_vectors':vectors}
        absolute=defaultdict(lambda:defaultdict(float));moments=defaultdict(lambda:defaultdict(float));totals=defaultdict(lambda:defaultdict(float))
        b.accumulate_vectors(vectors,absolute,moments,totals,{})
        stats=b.shape_statistics(result,absolute,moments,totals)
        first=(*channel,3)
        self.assertEqual(stats[first]['gen2'],0.)
        self.assertEqual(stats[first]['gen_gen'],0.)
        for row in stats.values():self.assertEqual(row['gen_gen'],row['gen2'])
        self.assertTrue(any(row['key']==list(first) and row['field']=='gen' and row['negative_numerator']<0
                            for row in result['shape_moment_roundoff']['cells']))
        _,direct,_=plot.normalized_shapes(result)
        self.assertAlmostEqual(direct[first]['gen2'],stats[first]['gen2'],delta=1e-30)
        # A materially inconsistent moment must still fail the existing guard.
        absolute[first]['gen2']=0.
        with self.assertRaisesRegex(ValueError,'beyond documented floating-point roundoff'):
            b.shape_statistics(result,absolute,moments,totals)

if __name__=='__main__':unittest.main()
