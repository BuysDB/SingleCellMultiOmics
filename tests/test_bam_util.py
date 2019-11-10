#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
import pysamiterators.iterators
from singlecellmultiomics.bamProcessing import sorted_bam_file
import os
"""
These tests check if the Molecule module is working correctly
"""

class TestSorted(unittest.TestCase):



    def test_write_to_sorted(self):
        write_path = './data/write_test.bam'
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            with sorted_bam_file(write_path, origin_bam=f) as out:
                for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
                    fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
                    fragment_class_args={'umi_hamming_distance':0},
                    pooling_method=0,
                    yield_invalid=True
                ):
                    molecule.write_pysam(out)

        self.assertTrue(os.path.exists(write_path))
        try:
            os.remove(write_path)
        except Exception as e:
            pass

    def test_write_to_read_grouped_sorted(self):
        write_path = './data/write_test_rg.bam'
        read_groups = set()
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            with sorted_bam_file(write_path, origin_bam=f,read_groups=read_groups) as out:
                for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
                    fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
                    fragment_class_args={'umi_hamming_distance':0},
                    pooling_method=0,
                    yield_invalid=True
                ):
                    molecule.write_pysam(out)
                    for frag in molecule:
                        read_groups.add( frag.get_read_group() )

        self.assertTrue(os.path.exists(write_path))
        try:
            os.remove(write_path)
        except Exception as e:
            pass


if __name__ == '__main__':
    unittest.main()
