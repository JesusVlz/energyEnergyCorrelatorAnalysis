import math
import unittest
from types import SimpleNamespace
from compare_calibrations import native_match, output_bin, weighted_profile
from veto_rule import reject_event
import numpy as np


class DiagnosticTests(unittest.TestCase):
    def test_output_boundaries(self):
        self.assertEqual([output_bin(x) for x in [39.99,40,79.99,80,120,160]],[0,1,1,2,3,4])

    def test_match_preserves_native_order_and_tolerances(self):
        jets=SimpleNamespace(refpt=[50.],refeta=[.1],refphi=[.2],ngen=2,
                             genpt=[50.5,50.],geneta=[.1,.1],genphi=[.2,.2])
        self.assertEqual(native_match(jets,0),0)
        jets.refpt=[-999.]
        self.assertIsNone(native_match(jets,0))

    def test_profile_exposes_empty_and_signed_weight_support(self):
        x=np.array([1.,1.]);y=np.array([1.,3.]);w=np.array([1.,1.])
        r=weighted_profile(x,y,w,np.array([0.,2.,4.]))
        self.assertEqual((r[0]['mean'],r[0]['width'],r[0]['neff']),(2.,1.,2.))
        self.assertIsNone(r[1]['mean'])
        self.assertIsNone(weighted_profile(x,y,np.array([-1.,-1.]),np.array([0.,2.]))[0]['mean'])

    def test_event_veto_exact_boundaries_and_nonzero(self):
        def jet(pt=16,cef=.3,nef=.4):
            return dict(pt=pt,eta=.2,phi=.3,charged_em=cef,neutral_em=nef)
        kwargs=dict(pt_key='pt',eta_key='eta',phi_key='phi')
        self.assertFalse(reject_event([jet(15)],lambda e,p:1,lambda j:True,**kwargs))
        self.assertTrue(reject_event([jet(15.000001)],lambda e,p:1,lambda j:True,**kwargs))
        self.assertFalse(reject_event([jet(16,.4,.5)],lambda e,p:1,lambda j:True,**kwargs))
        self.assertFalse(reject_event([jet()],lambda e,p:1,lambda j:False,**kwargs))
        self.assertTrue(reject_event([jet(10),jet()],lambda e,p:-1,lambda j:True,**kwargs))
        self.assertFalse(reject_event([jet()],lambda e,p:0,lambda j:True,**kwargs))

    def test_veto_no_invalid_map_fallback(self):
        jet=dict(pt=20,eta=.2,phi=.3,charged_em=.2,neutral_em=.3)
        with self.assertRaises(ValueError):
            reject_event([jet],lambda e,p:math.nan,lambda j:True,pt_key='pt',eta_key='eta',phi_key='phi')


if __name__=='__main__':unittest.main()
