import unittest
from joint_statistics import reduce_reports


def report(namespace, radius, weight):
    identity = [namespace, 1, 2, 3]
    return {"event_folds": [{"id": identity, "fold": 1}],
            "holdout_event_vectors": [{"id": identity, "fold": 1,
              "channels": [{"key": [namespace, radius, 0, 1, 2],
                            "gen": weight, "raw": weight, "single": weight,
                            "pair_supported": weight, "masked_single": 0}]}]}


def entry(result, radius1, radius2):
    first = result["columns"].index(["OO", radius1, 0, 1, 2, "raw"])
    second = result["columns"].index(["OO", radius2, 0, 1, 2, "raw"])
    return {(i,j): v for i,j,v in result["covariance_upper"]}.get(tuple(sorted((first, second))), 0)


class SharedEvents(unittest.TestCase):
    def test_cross_radius_shared_event(self):
        result = reduce_reports([report("OO", .4, 2), report("OO", .8, 3)])
        self.assertEqual(result["physical_events"], 1)
        self.assertEqual(entry(result, .4, .8), 6)
        self.assertEqual(entry(result, .4, .4), 4)

    def test_signed_event_weights_keep_correlation(self):
        result = reduce_reports([report("OO", .4, -2), report("OO", .8, 3)])
        self.assertEqual(entry(result, .4, .8), -6)

    def test_duplicate_input_rejected(self):
        with self.assertRaises(ValueError):
            reduce_reports([report("OO", .4, 2), report("OO", .4, 2)])

    def test_cross_radius_training_leakage_rejected(self):
        second = report("OO", .8, 3)
        second["event_folds"][0]["fold"] = 0
        with self.assertRaises(ValueError):
            reduce_reports([report("OO", .4, 2), second])

    def test_empty_event_retains_holdout_exposure(self):
        source = report("OO", .4, 2)
        source["holdout_event_vectors"][0]["channels"] = []
        result = reduce_reports([source])
        self.assertEqual(result["physical_events"], 1)
        self.assertEqual(result["contributing_physical_events"], 0)
        self.assertEqual(result["covariance_upper"], [])

    def test_independent_systems(self):
        result = reduce_reports([report("OO", .4, 2), report("pp", .4, 3)])
        self.assertEqual(result["physical_events"], 2)
        columns = result["columns"]
        for i,j,value in result["covariance_upper"]:
            self.assertEqual(columns[i][0], columns[j][0])


if __name__ == "__main__":
    unittest.main()
