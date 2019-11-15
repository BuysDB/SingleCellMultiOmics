#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
import pysamiterators.iterators
from singlecellmultiomics.bamProcessing import sorted_bam_file,write_program_tag
import os
import sys
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

            input_header = f.header.as_dict()
            write_program_tag(input_header,
                program_name='test_bam_util_test1',
                command_line = " ".join(sys.argv),
                version = singlecellmultiomics.__version__,
                description = f'a description'
                )

            write_program_tag(input_header,
                program_name='test_bam_util_test2',
                command_line = " ".join(sys.argv),
                version = singlecellmultiomics.__version__,
                description = f'a description'
                )
            #print([x for x in input_header['PG'] if not 'bwa mem' in x.get('CL','')])
            with sorted_bam_file(write_path, header=input_header,read_groups=read_groups) as out:
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

        # Now test if the program tag is there...
        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'test_bam_util_test1' in x.get('PN','')]) )
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'test_bam_util_test2' in x.get('PN','')]) )
            i =0

            # Test if the file has reads.
            for read in f:
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 293)
        try:
            os.remove(write_path)
        except Exception as e:
            pass


if __name__ == '__main__':
    unittest.main()
