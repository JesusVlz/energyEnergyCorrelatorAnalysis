import unittest
from pilot_input_audit import bounded_ranges, event_key, physical_id, partition, background_key, POOL_NAMESPACE


class PilotInputAuditTests(unittest.TestCase):
    def test_global_cap_does_not_expand_per_file(self):
        self.assertEqual(bounded_ranges([30097, 30000, 30000, 30000, 30000], 100000), [30097,30000,30000,9903,0])
        self.assertEqual(bounded_ranges([3,4], 20), [3,4])
        with self.assertRaises(ValueError): bounded_ranges([3], 0)

    def test_exact_64bit_identity_and_partition(self):
        ids = [[1,2,(1 << 53)+i] for i in range(100)]
        self.assertNotEqual(event_key(POOL_NAMESPACE, ids[0]),event_key(POOL_NAMESPACE, ids[1]))
        a = {tuple(k) for k in ids if partition(k,20261030)=='A'}
        b = {tuple(k) for k in ids if partition(k,20261030)=='B'}
        self.assertFalse(a & b)
        self.assertEqual(a | b, set(map(tuple,ids)))
        self.assertEqual(a, {tuple(k) for k in reversed(ids) if partition(k,20261030)=='A'})
        self.assertEqual(event_key('toy',[1,2,3]), 8856647007772020293)

    def test_same_background_stays_in_same_partition_across_physical_ids(self):
        signature = '1003:80daaf27abfea734:4e5d145bb4ac8308'
        self.assertEqual(partition([1,14,13542],20261030,signature), partition([1,53,50845],20261030,signature))
        self.assertIsNotNone(background_key(signature))
        empty = '0:0000000000000000:0000000000000000'
        self.assertIsNone(background_key(empty))
        self.assertEqual(partition([1,2,3],20261030,empty), partition([1,2,3],20261030))

    def test_invalid_ids_do_not_fall_back_to_entry(self):
        for identity in ([1,2], [1,2,-1], [1,2,1 << 64], [1,2,1.0], [True,2,3]):
            with self.assertRaises(ValueError): physical_id(identity)
        with self.assertRaises(ValueError): event_key('',[1,2,3])
        with self.assertRaises(ValueError): partition([1,2,3], -1)


if __name__ == '__main__': unittest.main()
