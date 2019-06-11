#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import itertools

import singlecellmultiomics.universalBamTagger.universalBamTagger as ut

"""
These tests check if the tagger is working correctly
"""

class TestMoleculeIteration(unittest.TestCase):

    def test_iteration(self):
        # Todo
        #ut.MoleculeIterator(a, umi_hamming_distance=1)
        pass

    def test_umi_hamming_iteration(self):
        pass

if __name__ == '__main__':
    unittest.main()
