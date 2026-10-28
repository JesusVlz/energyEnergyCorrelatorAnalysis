import copy
import unittest
from channel_contract import channels, validate_ratio


class RatioGate(unittest.TestCase):
    def setUp(self):
        self.oo = dict(system="OO", energy_tev=5.36, radius=0.4, axis="WTA",
                       jet_pt=[40, 80], track_pt_min=2, track_eta_max=2.4,
                       track_cut_version="synthetic_test_only", energy_power=1,
                       angular_edges=[0.01, 0.1, 0.4], normalization_range=[0.01, 0.4],
                       stage="final_normalized", input_hashes=["synthetic"],
                       calibration_hashes=["synthetic"], native_collection_verified=True,
                       correction_applicability_verified=True, independent_closure_verified=True,
                       covariance_verified=True, normalization_range_verified=True)
        self.pp = dict(self.oo, system="pp")

    def test_all_twelve_unique(self):
        self.assertEqual(len({x["channel"] for x in channels()}), 12)

    def test_matching_synthetic_metadata(self):
        self.assertTrue(validate_ratio(self.oo, self.pp))

    def test_incompatible_denominators(self):
        for key, value in (("radius", 0.8), ("energy_power", 2), ("jet_pt", [80, 120]),
                           ("normalization_range", [0.1, 0.4]), ("stage", "raw")):
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_ratio(self.oo, dict(self.pp, **{key: value}))

    def test_common_invalid_inputs_still_blocked(self):
        for key, value in (("correction_applicability_verified", False),
                           ("independent_closure_verified", False),
                           ("covariance_verified", False), ("input_hashes", []),
                           ("normalization_range", None), ("angular_edges", [0.4, 0.1])):
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_ratio(dict(self.oo, **{key: value}), dict(self.pp, **{key: value}))


if __name__ == "__main__":
    unittest.main()
