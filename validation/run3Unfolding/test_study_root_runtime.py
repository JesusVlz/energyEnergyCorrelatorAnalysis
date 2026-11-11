"""Real RooUnfold regression; require explicitly on the analysis host."""
import os
import unittest
import numpy as np
from study_root import error_enum,t_matrix,unfolded_result

class RootRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import ROOT
            ROOT.gROOT.SetBatch(True)
            if ROOT.gSystem.Load(os.environ.get('ROOUNFOLD_LIBRARY','libRooUnfold'))<0:
                raise RuntimeError('RooUnfold unavailable')
            cls.root=ROOT
        except (ImportError,RuntimeError) as error:
            if os.environ.get('EEC_REQUIRE_ROOT_TESTS')=='1': raise
            raise unittest.SkipTest(str(error))
    def test_rectangular_covariance(self):
        r=self.root
        m=r.TH1D('rt_m','',3,0,3);t=r.TH1D('rt_t','',2,0,2)
        c=r.TH2D('rt_c','',3,0,3,2,0,2)
        counts=np.array([[60.,10.],[20.,50.],[0.,10.]])
        for j in range(2): t.SetBinContent(j+1,100.)
        for i in range(3):
            m.SetBinContent(i+1,float(counts[i].sum()+5.))
            for j in range(2): c.SetBinContent(i+1,j+1,float(counts[i,j]))
        response=r.RooUnfoldResponse(m,t,c)
        solver=r.RooUnfoldBayes(response,m,3)
        solver.SetMeasuredCov(t_matrix(r,np.diag([80.,80.,20.])))
        h,v,api=unfolded_result(solver,error_enum(r))
        self.assertEqual(h.GetNbinsX(),2);self.assertEqual(v.shape,(2,2))
        self.assertTrue(np.all(np.isfinite(v)));self.assertNotEqual(float(v[0,1]),0.)
        self.assertIn(api['covariance_getter'],('Eunfold','Ereco'))

if __name__=='__main__': unittest.main()
