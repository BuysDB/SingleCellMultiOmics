#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from singlecellmultiomics.methylation import MethylationCountMatrix

"""
These tests check if the methylation module is working correctly
"""

class TestMethylationCountMatrix(unittest.TestCase):

    def test_matrix(self):
        mA = MethylationCountMatrix()
        mA['sample_A', ('chr1', 10, 20)][0] += 1
        mA['sample_A', ('chr2', 10, 20)][0] += 2

        self.assertTrue( len(mA.counts) == 1 )
        self.assertEqual(mA.counts['sample_A'].get(('chr1', 10, 20))[0], 1)
        self.assertEqual(mA.counts['sample_A'].get(('chr1', 10, 20))[1], 0)
        self.assertEqual(mA.counts['sample_A'].get(('chr2', 10, 20))[0], 2)
        self.assertEqual(mA.counts['sample_A'].get(('chr2', 10, 20))[1], 0)

        self.assertEqual(mA.get_bulk_frame('pd').loc[('chr1', 10, 20)].unmethylated , 1.0)
        self.assertEqual(mA.get_bulk_frame('pd').loc[('chr1', 10, 20)].methylated , 0)
        self.assertEqual(mA.get_bulk_frame('pd').loc[('chr1', 10, 20)].beta , 0 )

        mB = MethylationCountMatrix()
        mB['sample_A', ('chr2', 10, 20)][1] += 2
        mB['sample_B', ('chr2', 10, 20)][1] += 1
        mB['sample_B', ('chr2', 10, 20)][0] = 2

        mA.update(mB)
        self.assertEqual(len(mA.counts), 2)
        self.assertEqual(mA.get_frame('beta')[('chr2', 10, 20)]['sample_A'], 1)
        self.assertEqual(mA.get_frame('beta')[('chr2', 10, 20)]['sample_B'], 1 / 3)

        self.assertEqual(mA.get_bulk_frame()['n_samples'].max(), 2)
        self.assertEqual(mA.get_bulk_frame()['n_samples'].min(), 1)
        self.assertEqual(mA.get_bulk_frame()['methylated'].sum(), 3)
        self.assertEqual(mA.get_bulk_frame()['unmethylated'].sum(), 3)



if __name__ == '__main__':
    unittest.main()
