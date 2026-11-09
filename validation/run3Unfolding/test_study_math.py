import gzip
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from study_math import *
from study_root import CONTRACT_KEYS,assert_compatible,sha256,unfolded_result
from run3_study import validate_channel,verify_split


class MathTests(unittest.TestCase):
    def test_serialized_edges(self): self.assertEqual(physical_index([.7,1,2.4],float(np.float32(2.4))),2)
    def test_missing_edge(self):
        with self.assertRaises(ValueError): physical_index([40,80,120],39.9999)
    def test_no_subbin_interpolation(self):
        with self.assertRaises(ValueError): interval_index([40,80,120],[40,120])
    def test_invalid_edges(self):
        for e in ([1,1,2],[2,1],[0,np.nan]):
            with self.assertRaises(ValueError): edges(e)
    def test_serialization_degeneracy(self):
        with self.assertRaises(ValueError): physical_index([1,1+1e-10,2],1)
    def test_flat_layout(self):
        n,dr=flattened_layout([0,.1,.4,1,1.1,1.4,2],[5,40,80])
        self.assertEqual(n,3);np.testing.assert_allclose(dr,[0,.1,.4,1])
    def test_bad_blocks(self):
        with self.assertRaises(ValueError): flattened_layout([0,.1,1,1.2,2],[5,40,80])
    def test_bad_flat_dimensions(self):
        with self.assertRaises(ValueError): flattened_layout([0,.1,1,2],[5,40,80])
    def test_efficiency_and_fakes(self):
        c=np.array([[60.,10.],[20.,50.],[0.,10.]])
        t=np.array([100.,100.]);f=np.array([5.,10.,3.]);m=c.sum(axis=1)+f
        k,miss,fake=response_operator(c,m,t)
        np.testing.assert_allclose(miss,[20,30]);np.testing.assert_allclose(fake,f)
        np.testing.assert_allclose(forward_fold(k,t,fake,1),m)
    def test_not_unit_columns(self): self.assertEqual(response_operator([[40]],[40],[100])[0][0,0],.4)
    def test_bad_bookkeeping(self):
        with self.assertRaises(ValueError): response_operator([[10]],[5],[10])
    def test_roundoff_not_clipped(self):
        k,m,f=response_operator([[100.0001]],[100],[100]);self.assertLess(m[0],0);self.assertLess(f[0],0)
    def test_negative_response(self):
        with self.assertRaises(ValueError): response_operator([[-1]],[1],[1])
    def test_zero_truth(self):
        k,m,f=response_operator([[0]],[0],[0]);np.testing.assert_equal(k,[[0]])
    def test_rectangular_nullity(self):
        d=diagnostics(np.array([[1.,0,0],[0,1,0]]));self.assertEqual(d['truth_nullity'],1);self.assertEqual(d['condition_number'],'infinite')
    def test_identity_condition(self): self.assertEqual(diagnostics(np.eye(2))['condition_number'],1.)
    def test_cov_dimension(self):
        with self.assertRaises(ValueError): covariance(np.eye(3),2)
    def test_cov_symmetry(self):
        with self.assertRaises(ValueError): covariance([[1,.2],[0,1]])
    def test_indefinite_covariance(self):
        with self.assertRaises(ValueError): covariance([[1,2],[2,1]],check_psd=True)
    def test_singular_covariance(self): np.testing.assert_allclose(covariance([[1,1],[1,1]],check_psd=True),[[1,1],[1,1]])
    def test_refold_dimensions(self):
        with self.assertRaises(ValueError): forward_fold(np.eye(2),[1],[0,0],1)
    def test_refold_scale(self):
        with self.assertRaises(ValueError): forward_fold(np.eye(2),[1,1],[0,0],-1)
    def test_normalization_jacobian(self):
        y=np.array([2.,3.,5.]);w=np.array([.1,.2,.3]);m=np.ones(3,bool)
        f,v,j=normalize(y,np.eye(3),w,m);self.assertAlmostEqual(float(f@w),1.)
        numeric=np.empty_like(j)
        for k in range(3):
            dy=np.zeros(3);dy[k]=1e-5
            numeric[:,k]=(normalize(y+dy,np.eye(3),w,m)[0]-normalize(y-dy,np.eye(3),w,m)[0])/(2e-5)
        np.testing.assert_allclose(j,numeric,rtol=1e-8,atol=1e-9)
        self.assertAlmostEqual(float(w@v@w),0.)
    def test_nonpositive_norm(self):
        with self.assertRaises(ValueError): normalize([0,0],np.eye(2),[1,1],[1,1])
    def test_native_endpoints(self):
        m,e=normalization_bins([0,.005,.01,.1,.4,.8],[.008,.39])
        self.assertEqual(m.tolist(),[False,True,True,True,False]);self.assertEqual(e,[.005,.4])
    def test_ratio_joint_covariance(self):
        y=np.array([2.,4.]);v=np.array([[1,.3],[.3,2.]])
        r,vr=ratio(y,v,y,v,v);np.testing.assert_allclose(r,1);np.testing.assert_allclose(vr,0,atol=1e-15)
    def test_zero_denominator(self):
        with self.assertRaises(ValueError): ratio([1],[[1]],[0],[[1]])
    def test_one_candidate_per_matrix(self):
        jets=np.array([40.,80.,120.,160.]);dr=np.array([0.,.01,.1,.4,.8]);data={}
        truth=np.tile([1.,2.,3.,1.],3)
        for s in (0,1):
            data[f'truth_{s}']=truth;data[f'truth_var_{s}']=truth
            for k in range(1,21):
                data[f'u_{s}_{k}']=truth+np.tile([0,.02*(k-6),0,0],3)
                data[f'cov_{s}_{k}']=np.eye(12)
        rows,candidate=summarize_scan(data,jets,dr,[.008,.39])
        self.assertEqual(len(rows),120);self.assertEqual(candidate['iterations'],6)
        self.assertFalse(candidate['frozen']);self.assertEqual(data['chi2'].shape,(2,20))


class MatrixFake:
    def __init__(self,a): self.a=a
    def GetNrows(self): return len(self.a)
    def GetNcols(self): return len(self.a[0])
    def __getitem__(self,i): return self.a[i]

class HistogramFake:
    def GetNbinsX(self): return 2

class SolverNew:
    def __init__(self): self.called=[]
    def Hunfold(self,mode): self.called.append('Hunfold');return HistogramFake()
    def Eunfold(self,mode): self.called.append('Eunfold');return MatrixFake([[4.,1.],[1.,9.]])
    def GetMeasuredCov(self): raise AssertionError('Never use input covariance as result')

class SolverOld:
    def Hreco(self,mode): return HistogramFake()
    def Ereco(self,mode): return MatrixFake([[4.,1.],[1.,9.]])

class ContractTests(unittest.TestCase):
    def test_modern_result_api(self):
        s=SolverNew();h,v,api=unfolded_result(s,2)
        self.assertEqual(s.called,['Hunfold','Eunfold']);np.testing.assert_equal(v,[[4,1],[1,9]])
    def test_legacy_result_api(self): self.assertEqual(unfolded_result(SolverOld(),2)[2]['covariance_getter'],'Ereco')
    def test_no_measured_fallback(self):
        with self.assertRaises(RuntimeError): unfolded_result(object(),2)
    def test_output_covariance_space(self):
        s=SolverNew();s.Eunfold=lambda _:MatrixFake(np.eye(3))
        with self.assertRaises(ValueError): unfolded_result(s,2)
    def test_bad_power(self):
        with self.assertRaises(ValueError): validate_channel(dict(id='OO',system='OO',power=3,activity=[0,30]))
    def test_unsafe_name(self):
        with self.assertRaises(ValueError): validate_channel(dict(id='../OO',system='OO',power=1,activity=[0,30]))
    def test_partial_reference(self):
        with self.assertRaises(ValueError): validate_channel(dict(id='OO',system='OO',power=1,activity=[0,30],train_gen='x'))
    def test_truth_not_substituted(self):
        with self.assertRaises(ValueError): validate_channel(dict(id='OO',system='OO',power=1,activity=[0,30],truth_definition='ordinary_GenGen'))
    def test_partial_fake_exposure(self):
        with self.assertRaises(ValueError): validate_channel(dict(id='OO',system='OO',power=1,activity=[0,30],fake_scales={'A_to_B':1}))
    def test_wrong_system(self):
        a={k:[1] for k in CONTRACT_KEYS};b=dict(a);b['DataType']=[12]
        with self.assertRaises(ValueError): assert_compatible(a,b)
    def ledger(self,duplicate=False):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);ledger=p/'events.csv.gz'
            with gzip.open(ledger,'wt') as f: f.write('event_key,fold\na,0\n'+('a' if duplicate else 'b')+',1\n')
            hashes={'train':'abc','validation':'def'}
            proof=dict(product_sha256=hashes,ledgers=[dict(path=ledger.name,sha256=sha256(ledger))])
            (p/'proof.json').write_text(json.dumps(proof))
            return verify_split({'split_evidence':'proof.json'},p,hashes)
    def test_partition(self): self.assertEqual(self.ledger()['events'],[1,1])
    def test_duplicate_identity(self):
        with self.assertRaises(ValueError): self.ledger(True)


if __name__=='__main__': unittest.main()
