"""Synthetic rendering tests only; no physics/MC conclusions."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from study_math import summarize_scan
from draw_run3_study import draw

class PlotTests(unittest.TestCase):
    def make(self,p,scan):
        jets=np.array([5.,40.,80.,120.,160.,200.]);dr=np.array([0.,.01,.1,.4,.8]);n=20
        data=dict(jets=jets,dr=dr,response_0=np.eye(n),response_1=np.eye(n)*2)
        record=dict(synthetic=True,status='SCAN_EXECUTED_DIAGNOSTIC_ONLY' if scan else 'PREFLIGHT_ONLY',
            channel=dict(system='OO',activity=[0,30],power=1),truth_definition='response_marginal')
        if scan:
            truth=np.tile([1.,2.,3.,1.],5)
            for s in (0,1):
                data[f'truth_{s}']=truth;data[f'truth_var_{s}']=truth
                for k in range(1,21):
                    data[f'u_{s}_{k}']=truth+np.tile([0,.02*(k-6),0,0],5)
                    data[f'cov_{s}_{k}']=np.eye(n)
            _,record['candidate']=summarize_scan(data,jets,dr,[.008,.39])
        np.savez_compressed(p/'plot_data.npz',**data)
        (p/'manifest.json').write_text(json.dumps(record))
    def test_preflight_only_responses(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);self.make(p,False);draw(p,formats=('png',))
            self.assertEqual(len(list((p/'figures').glob('*.png'))),2)
    def test_all_plot_groups(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);self.make(p,True);draw(p,formats=('png',))
            self.assertEqual(len(list((p/'figures').glob('*.png'))),5)
    def test_no_placeholders(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);self.make(p,False);(p/'manifest.json').write_text(json.dumps({'status':'BLOCKED'}))
            with self.assertRaises(ValueError): draw(p,formats=('png',))
            self.assertFalse((p/'figures').exists())

if __name__=='__main__': unittest.main()
