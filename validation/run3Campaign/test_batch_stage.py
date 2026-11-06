import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('stage',Path(__file__).with_name('batch_stage.py'))
stage=importlib.util.module_from_spec(spec);spec.loader.exec_module(stage)

class PublicationTests(unittest.TestCase):
    def test_completion_is_published_last_and_binds_files(self):
        with tempfile.TemporaryDirectory() as work:
            directory=Path(work);(directory/'result.root').write_bytes(b'root fixture')
            commands=[]
            with patch.object(stage,'run',side_effect=lambda *args:commands.append(args)):
                stage.publish(directory,'root://example//output')
            copies=[c for c in commands if c[0]=='xrdcp']
            self.assertTrue(str(copies[-1][-1]).endswith('/stage_complete.json'))
            mkdir_commands=[c for c in commands if c[0]=='xrdfs']
            self.assertTrue(mkdir_commands)
            self.assertTrue(all(not str(c[-1]).endswith('/.') for c in mkdir_commands))
            self.assertTrue(all(not str(c[-1]).startswith('//') for c in mkdir_commands))
            receipt=json.loads((directory/'stage_complete.json').read_text())
            self.assertEqual(receipt['files']['result.root']['sha256'],stage.digest(directory/'result.root'))
            self.assertNotIn('stage_complete.json',receipt['files'])

    def test_upstream_changed_artifact_is_rejected(self):
        with tempfile.TemporaryDirectory() as work:
            target=Path(work)
            def fake_get(url,path):
                path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
                if url.endswith('stage_complete.json'):
                    path.write_text(json.dumps({'status':'COMPLETE','files':{'a.root':{'bytes':3,'sha256':'invalid'}}}))
                else:path.write_bytes(b'bad')
                return path
            with patch.object(stage,'get',side_effect=fake_get):
                with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                    stage.fetch_stage('root://example//output',target)

    def test_incomplete_upstream_is_rejected(self):
        with tempfile.TemporaryDirectory() as work:
            path=Path(work)/'stage_complete.json';path.write_text('{"status":"RUNNING"}')
            with patch.object(stage,'get',return_value=path):
                with self.assertRaisesRegex(ValueError,'incomplete'):
                    stage.fetch_stage('root://example//output',Path(work))

if __name__=='__main__':unittest.main()
