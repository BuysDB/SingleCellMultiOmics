#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
import pysamiterators.iterators
"""
These tests check if the Molecule module is working correctly
"""

class TestFragment(unittest.TestCase):

    def test_init(self):
        """Test if the fragment can be initialised"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for i,(R1,R2) in enumerate(pysamiterators.iterators.MatePairIterator(f)):
                frag = singlecellmultiomics.fragment.Fragment([R1,R2])
                self.assertIsNotNone(frag)

    def test_get(self):
        """Test if the fragment can be initialised"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for i,(R1,R2) in enumerate(pysamiterators.iterators.MatePairIterator(f)):
                frag = singlecellmultiomics.fragment.Fragment([R1,R2])
                if R1 is not None:
                    self.assertEqual(frag.get_R1(), R1 )
                    self.assertEqual(frag[0], R1 )
                if R2 is not None:
                    self.assertEqual(frag.get_R2(), R2 )
                    self.assertEqual(frag[1], R2 )

    def test_get_sample(self):
        """Test if the sample name of the fragment can be obtained"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for i,(R1,R2) in enumerate(pysamiterators.iterators.MatePairIterator(f)):
                frag = singlecellmultiomics.fragment.Fragment([R1,R2])
                self.assertTrue( frag.get_sample().startswith('A') )

    def test_set_sample(self):
        """Test if the sample name of the fragment can be changed"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for i,(R1,R2) in enumerate(pysamiterators.iterators.MatePairIterator(f)):
                frag = singlecellmultiomics.fragment.Fragment([R1,R2])
                frag.set_sample(f'TEST{i}')
                self.assertEqual(frag.get_sample(), f'TEST{i}' )

    def test_strand(self):
        """Test if the strand can be obtained (doesn't check validity)"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for i,(R1,R2) in enumerate(pysamiterators.iterators.MatePairIterator(f)):
                frag = singlecellmultiomics.fragment.Fragment([R1,R2])
                self.assertIn( frag.get_strand(), [None,0,1])
                if frag.is_mapped:
                    self.assertIn( frag.get_strand(), [0,1])


if __name__ == '__main__':
    unittest.main()
