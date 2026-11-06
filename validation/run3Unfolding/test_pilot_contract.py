import copy
import unittest
from pathlib import Path
from pilot_contract import edge_indices,float32,read_card,validate_card,validate_product_metadata

class PhysicalEdgeTests(unittest.TestCase):
    def test_serialization_only(self):
        self.assertEqual(edge_indices([.7,2.4],[float32(.7),1.,float32(2.4)]),[0,2])
    def test_missing_or_subdivided_edge(self):
        for requested in ([39.9999,80],[40,60,80],[10,40],[40,5020]):
            with self.assertRaises(ValueError): edge_indices(requested,[40,80,120,160])
    def test_invalid_edges(self):
        for available in ([40,40,80],[80,40],[40,float('nan')]):
            with self.assertRaises(ValueError): edge_indices([40],available)
    def test_artifact_cards(self):
        import json
        folder=Path(__file__).parent/'fixtures/pilot_cards'
        manifest=json.loads((folder/'manifest.json').read_text())
        self.assertEqual(len(manifest['tasks']),12)
        for task in manifest['tasks']:
            card=read_card(folder/task['card'])
            self.assertTrue(validate_card(card,task['system'],task['power'],task['activity_edges']))
            self.assertTrue(validate_product_metadata(card,copy.deepcopy(card),task['power']))
            for key,value in [('MinTrackPtCut',[.5]),('Run3JetID',[0]),('JetAxis',[0]),('MinJetPtUnfoldingReco',[80]),('WeightExponent',[3]),('CentralityBinEdges',[-1,90])]:
                actual=copy.deepcopy(card);actual[key]=value
                with self.assertRaises(ValueError): validate_product_metadata(card,actual,task['power'])
            actual=copy.deepcopy(card);del actual['Run3JetID']
            with self.assertRaises(ValueError): validate_product_metadata(card,actual,task['power'])
            actual=copy.deepcopy(card);actual['JetPtBinEdgesEEC']=[40,80,120,160,200,300,500,5020]
            with self.assertRaises(ValueError): validate_product_metadata(card,actual,task['power'])

if __name__=='__main__': unittest.main()
