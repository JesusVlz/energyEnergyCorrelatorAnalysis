import argparse
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('closure_worker',Path(__file__).with_name('closure_worker.py'))
worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
REPO=Path(__file__).resolve().parents[1]


class FrozenWorkerIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import ROOT
        ROOT.gROOT.SetBatch(True);cls.ROOT=ROOT

    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory(prefix='oo_closure_worker_');self.addCleanup(self.temporary.cleanup)
        self.base=Path(self.temporary.name);self.analysis=self.base/'analysis';self.analysis.mkdir()
        (self.analysis/'cards').mkdir();(self.analysis/'mixing_pool').mkdir();(self.analysis/'mixingFileList').mkdir();(self.analysis/'condor').mkdir()
        shutil.copy2(REPO/'condor/checkOOMCOutput.C',self.analysis/'condor/checkOOMCOutput.C')
        (self.analysis/'eecAnalysis').write_text('native executable fixture; execution is explicitly mocked')
        (self.analysis/'campaign_provenance.json').write_text(json.dumps({'source_commit':'testfixture','executable_sha256':worker.sha256(self.analysis/'eecAnalysis')}))
        self.mapping={};self.publications={};self.commands=[];self.variant=None
        entries=[]
        for index in (1,2):
            path=self.base/('source_pool_'+str(index)+'.root');self.source_root(path,1)
            metadata=worker.root_info(path);url='root://test//pool_'+str(index)+'.root';self.mapping[url]=path
            entries.append(dict(metadata,original_path='/eos/test/pool_'+str(index)+'.root',input_url=url,local_name='pool_%04d.root'%index,use_entries=1))
        self.signal=self.base/'signal.root';self.source_root(self.signal,3);self.mapping['root://test//signal.root']=self.signal
        signal_metadata=dict(worker.root_info(self.signal),path='/eos/test/signal.root')
        headers=self.analysis/'mixing_pool/signal_headers.manifest.json';headers.write_text(json.dumps({'files':[signal_metadata]}))
        audit=self.analysis/'mixing_pool/source_audit.manifest.json';audit.write_text('{"source":"existing pure HIJING"}')
        pool_list=self.analysis/'mixingFileList/OO2025_HIJING_MC_pool_campaign.txt';pool_list.write_text('\n'.join(row['original_path'] for row in entries)+'\n')
        manifest={'schema_version':1,'global_event_cap':2,'seed':20261030,'partition_policy':1,
                  'holdout_policy':'fixed_reserved_donors_exclude_overlapping_signals','reserved_pool_scope':'all scanned source IDs and observed background fingerprints',
                  'signal_max_entries':3,'entries':entries,'frozen_pool_list':str(pool_list.relative_to(self.analysis)),
                  'pool_list_sha256':worker.sha256(pool_list),'source_audit_sha256':worker.sha256(audit),'signal_headers_sha256':worker.sha256(headers)}
        (self.analysis/'mixing_pool/manifest.json').write_text(json.dumps(manifest))
        settings={'OOClosurePilot':1,'OOClosurePartitionPolicy':1,'OOClosureMaxSignalEvents':3,'OOClosureMaxPoolEvents':2,'OOClosureSeed':20261030,
                  'OOClosureOutputPrefix':'pilot','DataType':10,'McCorrelationType':3,'BackgroundMethods':6,'EventMatchForMixing':2,
                  'BinnedMixing':0,'MegaSkimMixing':0,'SubeventCut':2,'JetAxis':1,'JetRadius':.4,'HFEnergyShift':0,'SmearDeltaR':0,
                  'SmearEnergyWeight':0,'DisableTrackPairEfficiencyCorrection':1,'WeightExponent':'1 2','MixingFileList':str(pool_list.relative_to(self.analysis)),
                  'MaxEventsPerFile':0}
        self.card=self.analysis/'cards/GenGen.input';self.card.write_text('\n'.join(key+' '+str(value) for key,value in settings.items())+'\n')
        self.args=argparse.Namespace(analysis_dir=str(self.analysis),work_dir=str(self.base/'work'),card='cards/GenGen.input',
            input_url='root://test//signal.root',output_url='root://test//outputs/sample.root',stem='sample',archive_sha256='a'*64,
            mode='GenGen',data_type=10,preflight_only=False)

    def source_root(self,path,count):
        ROOT=self.ROOT;file=ROOT.TFile(str(path),'RECREATE')
        from array import array
        value=array('f',[1])
        for directory,name in [('hiEvtAnalyzer','HiTree'),('HiGenParticleAna','hi'),('skimanalysis','HltTree'),('akCs4PFJetAnalyzer','t')]:
            file.cd();file.mkdir(directory).cd();tree=ROOT.TTree(name,name);tree.Branch('fixture',value,'fixture/F')
            for _ in range(count):tree.Fill()
            tree.Write()
        file.Close()

    def output_root(self,path,jet_entries=0):
        ROOT=self.ROOT;from array import array
        file=ROOT.TFile(str(path),'RECREATE');file.mkdir('JCard').cd()
        for key,value in [('DataType',10),('TriggerSelection',0),('McCorrelationType',3)]:
            vector=ROOT.TVectorT('float')(1);vector[0]=value;vector.Write(key)
        file.cd();events=ROOT.TH1D('nEvents','nEvents',12,-.5,11.5)
        for x in (0,1,8):events.Fill(x,3)
        events.Write()
        for name in ('vertexZ','pthat','vertexZweighted','centrality','centralityWeighted','pthatWeighted'):
            histogram=ROOT.TH1D(name,name,10,0,10);histogram.Fill(1);histogram.Write()
        for name in ('inclusiveJet','energyEnergyCorrelator'):
            histogram=ROOT.THnSparseD(name,name,1,array('i',[10]),array('d',[0]),array('d',[100]))
            if name=='inclusiveJet':
                for _ in range(jet_entries):histogram.Fill(array('d',[50]))
            histogram.Write()
        file.Close()

    def ledger(self,path,rows,fields):
        with path.open('w',newline='') as output:
            writer=csv.DictWriter(output,fieldnames=fields);writer.writeheader();writer.writerows(rows)

    def fake_analysis(self,command):
        card=worker.settings(Path(command[2]).read_text());prefix=Path(card['OOClosureOutputPrefix'][0]);self.output_root(Path(command[3]))
        matched=self.variant=='false_empty_root'
        donors=[dict(donor_key=str(k),run='1',lumi='1',event=str(i),eligible='1',pool=group,background_key='0',generator_weight='1',nh='1',vz='0',uses='1' if matched else '0')
                for k,i,group in [(2,1,'A'),(3,2,'B')]]
        self.ledger(Path(str(prefix)+'.pool.csv'),donors,list(donors[0]))
        partners=[]
        for i in (1,2,99):
            row=dict(signal_key=str(100+i),run='1',lumi='1',event=str(i),signal_background_key='0',event_weight='1',inspected='0' if i!=99 else '2',
                     status='failed',reason='signal_in_reserved_donor_physical_id' if i!=99 else 'no_pool_A',vz='0',nh='1' if matched else '100')
            for partner in (1,2):
                for field in ('key','run','lumi','event','entry','vz','nh','weight'):row['donor'+str(partner)+'_'+field]=''
            if self.variant=='reserved_leak' and i==1:row['status']='matched';row['reason']='matched'
            if self.variant=='partial_failure' and i==1:row['donor1_key']='2'
            if matched and i==99:
                row['status']='matched';row['reason']='matched'
                for partner,donor in enumerate(donors,1):
                    for field in ('key','run','lumi','event','vz','nh'):row['donor'+str(partner)+'_'+field]=donor['donor_key' if field=='key' else field]
                    row['donor'+str(partner)+'_entry']=str(partner-1);row['donor'+str(partner)+'_weight']='1'
            partners.append(row)
        self.ledger(Path(str(prefix)+'.partners.csv'),partners,list(partners[0]))
        jets=[{'signal_key':row['signal_key'],'jet_index':'0','status':row['status']} for row in partners]
        self.ledger(Path(str(prefix)+'.jets.csv'),jets,list(jets[0]))
        contributions=[]
        if self.variant=='failed_contribution':contributions=[dict(signal_key='101',jet_index='0',power='1',value='1')]
        self.ledger(Path(str(prefix)+'.contributions.csv'),contributions,['signal_key','jet_index','power','value'])

    def execute(self,command,cwd):
        self.commands.append(command)
        if command[0]=='xrdcp':
            source,destination=command[-2:]
            if source.startswith('root://'):shutil.copy2(self.mapping[source],destination)
            else:self.publications[destination]=Path(source).read_bytes()
            return 0
        if command[0].endswith('/eecAnalysis'):self.fake_analysis(command);return 0
        if command[0]=='root':return subprocess.run(command,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).returncode
        raise AssertionError(command)

    def test_actual_root_headers_checker_and_mocked_transfer_full_lifecycle(self):
        result=worker.run_job(self.args,execute=self.execute)
        self.assertTrue(result['empty_common_exposure']);self.assertEqual(result['ledger_validation']['reserved_signal_failures'],2)
        self.assertEqual(len(self.publications),6)
        self.assertEqual(set(self.publications),{'root://test//outputs/sample'+suffix for suffix in ('.root',*worker.SUFFIXES,'.run.json')})
        metadata=json.loads(self.publications['root://test//outputs/sample.run.json'])
        self.assertEqual(metadata['source']['source_entries'],3);self.assertEqual(metadata['partition_policy'],1)
        frozen=worker.settings(self.card.read_text());runtime=worker.settings((Path(self.args.work_dir)/'runtime.input').read_text())
        self.assertEqual({key for key in frozen if frozen[key]!=runtime[key]},{'MixingFileList','OOClosureOutputPrefix'})
        self.assertEqual(metadata['artifacts']['sample.root']['sha256'],hashlib_sha(self.publications['root://test//outputs/sample.root']))

    def test_legacy_macro_keeps_nonempty_common_exposure_guard(self):
        output=self.base/'empty_legacy.root';self.output_root(output)
        invocation='condor/checkOOMCOutput.C('+json.dumps(str(output))+',"GenGen",10)'
        result=subprocess.run(['root','-l','-b','-q',invocation],cwd=self.analysis,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('inclusiveJet has no entries',result.stdout)

    def test_wrong_staged_source_rejected_before_analysis_or_publication(self):
        self.mapping['root://test//pool_1.root']=self.signal
        with self.assertRaisesRegex(ValueError,'Staged ROOT size differs'):
            worker.run_job(self.args,execute=self.execute)
        self.assertFalse(self.publications)
        self.assertFalse(any(command[0].endswith('/eecAnalysis') for command in self.commands))

    def test_preflight_does_not_stage_or_execute(self):
        self.args.preflight_only=True;result=worker.run_job(self.args,execute=self.execute)
        self.assertEqual(result['status'],'PREFLIGHT_PASS');self.assertFalse(self.commands)

    def test_frozen_hash_and_full_signal_cap_fail_before_staging(self):
        self.card.write_text(self.card.read_text().replace('OOClosureMaxSignalEvents 3','OOClosureMaxSignalEvents 2'))
        with self.assertRaisesRegex(ValueError,'truncates'):worker.run_job(self.args,execute=self.execute)
        self.assertFalse(self.commands)

    def test_partial_failure_reserved_leak_failed_contribution_and_false_empty_root_rejected(self):
        for variant in ('partial_failure','reserved_leak','failed_contribution','false_empty_root'):
            with self.subTest(variant=variant):
                self.variant=variant;self.args.work_dir=str(self.base/variant);self.publications.clear()
                with self.assertRaises(ValueError):worker.run_job(self.args,execute=self.execute)
                self.assertFalse(self.publications)

    def test_frozen_pool_list_tampering_rejected(self):
        (self.analysis/'mixingFileList/OO2025_HIJING_MC_pool_campaign.txt').write_text('/eos/unapproved/root.root\n')
        with self.assertRaisesRegex(ValueError,'hash mismatch'):worker.run_job(self.args,execute=self.execute)
        self.assertFalse(self.commands)


def hashlib_sha(data):
    import hashlib
    return hashlib.sha256(data).hexdigest()


if __name__=='__main__':unittest.main()
