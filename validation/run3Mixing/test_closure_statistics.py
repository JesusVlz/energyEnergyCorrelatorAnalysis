import math
import unittest
from closure_statistics import (derive, comparison, joint_moments, weighted_pairs,
                                SAME, JM1, JM2, M11, M22, M12, PP, PH, HP, HH, cluster_vectors)


class ClosureStatistics(unittest.TestCase):
    def toy(self, power=1, weight=1., py=(1., 1., 1.), hy=(1., 1., 1., 1.),
            donor1=(1., 1., 1., 1.), donor2=(1., 1., 1., 1.)):
        raw = {(SAME, PP): weighted_pairs(py, power=power),
               (SAME, PH): weighted_pairs(py, hy, power),
               (SAME, HH): weighted_pairs(hy, power=power),
               (JM1, PH): weighted_pairs(py, donor1, power),
               (JM1, HH): weighted_pairs(hy, donor1, power),
               (JM2, HP): weighted_pairs(py, donor2, power),
               (JM2, HH): weighted_pairs(hy, donor2, power),
               (M11, HH): weighted_pairs(donor1, power=power),
               (M22, HH): weighted_pairs(donor2, power=power),
               (M12, HH): weighted_pairs(donor1, donor2, power)}
        return {key: weight * value for key, value in raw.items()}

    def test_unit_counting_three_signal_four_background(self):
        raw = self.toy()
        self.assertEqual(raw[JM1, PH] + raw[JM1, HH], 28)
        self.assertEqual(raw[M11, HH], 6)
        self.assertEqual(raw[M12, HH], 16)
        result = derive(raw)
        self.assertEqual(result['Bmix'], 18)
        self.assertEqual(result['Smix'], 3)
        self.assertEqual(result['truth_sum_residual'], 0)
        self.assertEqual(result['mix_PMavg'], result['PH'])
        self.assertEqual(result['mix_HMavg_minus_cross'], 0)
        self.assertEqual(result['mix_same_avg'], result['HH'])

    def test_weighted_pairs_n1_n2(self):
        for power in (1, 2):
            py, hy = (2., 3., 5.), (1., 2., 4., 6.)
            result = derive(self.toy(power, py=py, hy=hy, donor1=hy, donor2=hy))
            self.assertAlmostEqual(result['Smix'], weighted_pairs(py, power=power))
            self.assertAlmostEqual(result['Bmix'], result['Btrue'])
            self.assertEqual(result['truth_sum_residual'], 0)

    def test_projected_ph_counted_once(self):
        result = derive({(SAME, PP): 2., (SAME, PH): 3., (SAME, HP): 5., (SAME, HH): 7.})
        self.assertEqual(result['PH'], 8.)
        self.assertEqual(result['T'], 17.)
        self.assertEqual(result['Btrue'], 15.)
        self.assertEqual(result['truth_sum_residual'], 0.)

    def test_donor_exchange_and_event_weight(self):
        first = derive(self.toy(2, -3., donor1=(1., 2.), donor2=(3., 5., 7.)))
        second = derive(self.toy(2, -3., donor1=(3., 5., 7.), donor2=(1., 2.)))
        for name in ('Bmix', 'Smix', 'mix_PMavg', 'mix_HMavg_minus_cross', 'mix_same_avg'):
            self.assertAlmostEqual(first[name], second[name])
        unweighted = derive(self.toy(2, donor1=(1., 2.), donor2=(3., 5., 7.)))
        for name in first:
            self.assertAlmostEqual(first[name], -3. * unweighted[name])

    def test_empty_cones_valid_zero(self):
        result = derive(self.toy(py=(), hy=(), donor1=(), donor2=()))
        self.assertTrue(all(value == 0. for value in result.values()))

    def test_joint_residual_uses_covariance(self):
        vectors = [derive(self.toy(weight=2.)), derive(self.toy(weight=3.))]
        totals, cov = joint_moments(vectors)
        check = comparison(totals, cov, 'Bmix', 'Btrue')
        self.assertEqual(check['residual'], 0.)
        self.assertEqual(check['residual_error'], 0.)
        self.assertEqual(check['ratio'], 1.)
        self.assertEqual(check['ratio_error'], 0.)
        self.assertEqual(cov['Btrue', 'Btrue'], 18. ** 2 * (2. ** 2 + 3. ** 2))

    def test_repeated_signal_backgrounds_clustered(self):
        first, second = derive(self.toy(weight=2.)), derive(self.toy(weight=3.))
        clustered = cluster_vectors({'large64bit:18446744073709551600': first,
                                    'large64bit:18446744073709551601': second},
                                   {'large64bit:18446744073709551600': 'background:100',
                                    'large64bit:18446744073709551601': 'background:100'})
        total, cov = joint_moments(clustered.values())
        self.assertEqual(len(clustered), 1)
        self.assertEqual(cov['PP', 'PP'], total['PP'] ** 2)

    def test_invalid_ratios_null(self):
        totals, cov = joint_moments([derive({})])
        check = comparison(totals, cov, 'Bmix', 'Btrue')
        self.assertFalse(check['valid_ratio'])
        self.assertIsNone(check['ratio'])
        self.assertIsNone(check['ratio_error'])

    def test_poisson_events_not_pairs(self):
        vector = derive(self.toy(weight=2.))
        totals, cov = joint_moments([vector])
        self.assertEqual(cov['PP', 'PP'], totals['PP'] ** 2)
        self.assertNotEqual(cov['PP', 'PP'], totals['PP'])


if __name__ == '__main__':
    unittest.main()
