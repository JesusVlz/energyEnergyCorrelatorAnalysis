import argparse
import importlib.util
import json
from pathlib import Path
import tarfile
import subprocess
import tempfile
import unittest
from unittest.mock import patch

MODULE_PATH = Path(__file__).with_name('prepareOOMC.py')
spec = importlib.util.spec_from_file_location('prepareOOMC', MODULE_PATH)
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)
REPO = MODULE_PATH.resolve().parents[1]


class CampaignPayloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = Path(self.tmp.name) / 'repo'
        self.repo.mkdir()
        self.relative = Path('jetVetoMaps/OO2025/test.root')
        self.map = self.repo / self.relative
        self.map.parent.mkdir(parents=True)
        self.map.write_bytes(b'fixed ROOT map bytes for payload preservation')
        self.sidecar = self.map.with_suffix('.metadata.json')
        self.sidecar.write_text(json.dumps({'filename':self.map.name, 'histogram_key':'jetvetomap_all',
                                           'sha256':prepare.sha256(self.map)}))
        self.card = 'ApplyJetVeto 1\nJetVetoMapFile '+str(self.relative)+'\nJetVetoMapHistogram jetvetomap_all\n'

    def test_worker_archive_contains_enabled_map_and_pinned_sidecar(self):
        files = prepare.veto_map_payload_files(self.repo, {'GenGen':self.card})
        snapshot = Path(self.tmp.name) / 'snapshot'
        prepare.copy_veto_map_payload(self.repo,snapshot,files)
        archive = Path(self.tmp.name) / 'payload.tar.gz'
        with tarfile.open(archive,'w:gz') as stream:stream.add(snapshot,arcname='analysis')
        with tarfile.open(archive) as stream:
            self.assertEqual(stream.extractfile('analysis/'+str(self.relative)).read(),self.map.read_bytes())
            metadata=json.load(stream.extractfile('analysis/'+str(self.relative.with_suffix('.metadata.json'))))
        self.assertEqual(metadata['sha256'],prepare.sha256(self.map))

    def test_missing_referenced_map_fails_before_creating_campaign(self):
        base = (REPO/'cardEECOOMC.input').read_text()
        (self.repo/'cardEECOOMC.input').write_text(base)
        for name in ('Makefile','eecAnalysis.cxx','condor/runOOMCJob.sh','condor/checkOOMCOutput.C'):
            target=self.repo/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_text('fixture')
        for name in ('src','jetEnergyCorrections','trackCorrectionTables'):(self.repo/name).mkdir()
        args=argparse.Namespace(repo=str(self.repo),sample='OOMC',input_dir='/eos/mock/input',submit_dir='/afs/mock/nonexistent/campaign',
            output_dir='/eos/mock/output',cmssw_base=self.tmp.name,scram_arch='el9_amd64_gcc12',tag='missingmap',modes=['GenGen'],
            build_jobs=1,accounting_group=None,list_only=False,background_methods=2,match_eec_binning=False,base_card=None)
        with patch.object(prepare,'parse_args',return_value=args), patch.object(prepare,'discover_inputs',return_value=[{}]), \
             patch.object(prepare,'git_output',return_value='fixture'), patch.object(Path,'mkdir') as mkdir:
            with self.assertRaisesRegex(ValueError,'Missing or empty referenced veto map'):
                prepare.main()
            mkdir.assert_not_called()

    def test_sidecar_checksum_and_presence_are_required(self):
        self.map.write_bytes(b'changed map')
        with self.assertRaisesRegex(ValueError,'checksum differs'):prepare.veto_map_payload_files(self.repo,{'GenGen':self.card})
        self.sidecar.unlink()
        with self.assertRaisesRegex(ValueError,'metadata sidecar'):prepare.veto_map_payload_files(self.repo,{'GenGen':self.card})

    def test_disabled_veto_does_not_require_a_map(self):
        self.assertEqual(prepare.veto_map_payload_files(self.repo,{'GenGen':'ApplyJetVeto 0\n'}),[])

    def test_custom_card_and_perpendicular_option_preserve_production_settings(self):
        base=(REPO/'cardEECOOMC.input').read_text().replace('CentralityBinEdges  0 30 50 90','CentralityBinEdges -1 100')
        card=prepare.card_text(base,'GenGen',2,'OOMC')
        settings={line.split('#')[0].split()[0]:line.split('#')[0].split()[1:] for line in card.splitlines() if line.split('#')[0].strip()}
        self.assertEqual(settings['CentralityBinEdges'],['-1','100'])
        self.assertEqual(settings['McCorrelationType'],['3'])
        self.assertEqual(settings['BackgroundMethods'],['2'])
        self.assertEqual(settings['Run3JetID'],['1'])
        self.assertEqual(settings['ApplyJetVeto'],['1'])
        self.assertNotIn('OOClosurePilot',settings)
        self.assertNotIn('OOClosureMaxSignalEvents',settings)
        with patch('sys.argv',['prepareOOMC.py','--base-card','condor/cards/integrated.input','--background-methods','2']):
            args=prepare.parse_args()
        self.assertEqual(args.base_card,'condor/cards/integrated.input')
        self.assertEqual(args.background_methods,2)


class ReservedPoolCampaignTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.repo=Path(self.tmp.name)
        signal=self.repo/'signal.root';signal.write_bytes(b'signal fixture')
        self.inputs=[{'path':str(signal),'size_bytes':signal.stat().st_size}]
        self.signal_manifest=self.repo/'signal_headers.json'
        self.signal_manifest.write_text(json.dumps({'status':'PASS_REQUIRED_HEADERS_ONLY','files':[
            {'path':str(signal),'size_bytes':signal.stat().st_size,'source_entries':30000,'status':'PASS_REQUIRED_HEADERS'}]}))
        source=[];paths=[];remaining=100000
        for i,count in enumerate((30097,29291,29291,30000,9999),1):
            path=self.repo/('donor'+str(i)+'.root');path.write_bytes(('donor fixture '+str(i)).encode());paths.append(str(path))
            used=min(count,remaining);remaining-=used
            source.append({'path':str(path),'size_bytes':path.stat().st_size,'source_entries':count,'root_uuid':'fixture-'+str(i),
                           'reco_veto_required_branches_present':True,'audit':{'entry_range_scanned':[0,used]}})
        self.pool_list=self.repo/'pool.txt';self.pool_list.write_text('\n'.join(paths)+'\n')
        self.pool_manifest=self.repo/'pool.manifest.json'
        self.pool_manifest.write_text(json.dumps({'input_pool_list_sha256':prepare.sha256(self.pool_list),'pool_files':source}))
        self.args=argparse.Namespace(background_methods=6,pool_events=100000,mixing_pool_list=str(self.pool_list),
            mixing_pool_manifest=str(self.pool_manifest),signal_input_manifest=str(self.signal_manifest))
        text=(REPO/'validation/run3Mixing/cards/OO_GenGen_R04_closure.input').read_text()
        text=prepare.replace_card_scalar(text,'OOClosureMaxSignalEvents',30000)
        if 'OOClosurePartitionPolicy' in prepare.card_settings(text):text=prepare.replace_card_scalar(text,'OOClosurePartitionPolicy',1)
        else:text+='\nOOClosurePartitionPolicy 1\n'
        self.cards={'GenGen':prepare.card_text(text,'GenGen',6,'OOMC')}

    def plan(self,cards=None,sample='OOMC',modes=None):
        with patch.object(prepare,'eos_url',side_effect=lambda path:'root://fixture//'+path.name):
            return prepare.mixing_pool_plan(self.repo,cards or self.cards,sample,modes or ['GenGen'],self.inputs,self.args)

    def test_fixed_100k_prefix_freezes_only_four_files_and_both_audits(self):
        plan=self.plan();manifest=plan['manifest']
        self.assertEqual([row['use_entries'] for row in manifest['entries']],[30097,29291,29291,11321])
        self.assertEqual(manifest['partition_policy'],1)
        self.assertEqual(manifest['signal_max_entries'],30000)
        self.assertEqual(sum(row['use_entries'] for row in manifest['entries']),100000)
        self.assertEqual([row['local_name'] for row in manifest['entries']],['pool_0001.root','pool_0002.root','pool_0003.root','pool_0004.root'])
        destination=self.repo/'snapshot';prepare.freeze_mixing_pool(plan,destination)
        frozen=json.loads((destination/'mixing_pool/manifest.json').read_text())
        self.assertEqual(frozen,manifest)
        self.assertEqual(prepare.sha256(destination/'mixing_pool/source_audit.manifest.json'),prepare.sha256(self.pool_manifest))
        self.assertEqual(prepare.sha256(destination/'mixing_pool/signal_headers.manifest.json'),prepare.sha256(self.signal_manifest))
        self.assertEqual(prepare.sha256(destination/manifest['frozen_pool_list']),manifest['pool_list_sha256'])
        self.assertEqual(prepare.sha256(self.pool_list),manifest['source_pool_list_sha256'])
        paths=(destination/manifest['frozen_pool_list']).read_text().splitlines()
        self.assertEqual(len(paths),4)
        self.assertNotIn('donor5.root',' '.join(paths))

    def test_method6_rejects_other_systems_modes_and_nonreserved_policy(self):
        for sample,modes in [('ppMC',['GenGen']),('OOMC',['RecoReco']),('OOMC',['GenGen','Response'])]:
            with self.assertRaisesRegex(ValueError,'restricted'):self.plan(sample=sample,modes=modes)
        for key,value in [('OOClosurePartitionPolicy',0),('OOClosureSeed',12),('MinPairTrackPt',.5)]:
            cards={'GenGen':prepare.replace_card_scalar(self.cards['GenGen'],key,value)}
            with self.assertRaisesRegex(ValueError,'requires'):self.plan(cards=cards)

    def test_full_signal_scope_and_actual_audited_pool_bound_required(self):
        cards={'GenGen':prepare.replace_card_scalar(self.cards['GenGen'],'OOClosureMaxSignalEvents',10000)}
        with self.assertRaisesRegex(ValueError,'complete signal file'):self.plan(cards=cards)
        m=json.loads(self.pool_manifest.read_text());m['pool_files'][3]['audit']['entry_range_scanned']=[0,100]
        self.pool_manifest.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,'beyond the bounded'):self.plan()

    def test_changed_donor_list_or_signal_bytes_are_rejected(self):
        self.pool_list.write_text(self.pool_list.read_text()+'# changed list\n')
        with self.assertRaisesRegex(ValueError,'changed since'):self.plan()
        self.inputs[0]['size_bytes']+=1
        with self.assertRaisesRegex(ValueError,'Signal input is missing or changed'):self.plan()



class OptimizedSnapshotTests(unittest.TestCase):
    def test_controlled_o2_compiler_override_preserves_existing_flags_for_all_objects(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            (folder/'Makefile').write_text('CXX = g++\nCXXFLAGS = -DKEEP_EXISTING_FLAGS\nall: first.o second.o\n')
            for name in ('first','second'):(folder/(name+'.cc')).write_text('int '+name+'(){return 0;}\n')
            result=subprocess.run(['make','-n','CXX=g++ -O2'],cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
            commands=result.stdout.splitlines()
            self.assertEqual(len(commands),2)
            self.assertTrue(all(line.startswith('g++ -O2 ') and '-DKEEP_EXISTING_FLAGS' in line for line in commands))
            self.assertNotIn('fast-math',result.stdout)

    def test_legacy_build_has_no_override_and_opted_build_records_exact_command(self):
        for optimize in (False,True):
            with tempfile.TemporaryDirectory() as temporary:
                folder=Path(temporary);exe=folder/'eecAnalysis';exe.write_text('fixture');exe.chmod(0o755)
                with patch.object(prepare.subprocess,'run',return_value=subprocess.CompletedProcess([],0)) as run:
                    prepare.build_snapshot(folder,folder,'el9_amd64_gcc12',2,folder/'build.log',optimize=optimize)
                script=run.call_args[0][0][-1]
                log=(folder/'build.log').read_text()
                self.assertEqual("CXX=g++ -O2" in script,optimize)
                self.assertEqual("CXX=g++ -O2" in log,optimize)
                self.assertNotIn('fast-math',script)



if __name__=='__main__':unittest.main()
