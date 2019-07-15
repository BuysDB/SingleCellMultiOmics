#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
"""
These tests check if the Molecule module is working correctly
"""

class TestMolecule(unittest.TestCase):

    def test_consensus(self):
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.Fragment)

            for molecule in iter(it):
                #print(molecule.get_sample())
                if  molecule.get_sample()=='APKS3-P19-1-1_91':
                    break

            self.assertEqual(''.join( list(molecule.get_consensus().values()) ), 'AGTTAGATATGGACTCTTCTTCAGACACTTTGTTTAAATTTTAAATTTTTTTCTGATTGCATATTACTAAAAATGTGTTATGAATATTTTCCATATCATTAAACATTCTTCTCAAGCATAACTTTAAATAACTATAGAAAATTTACGCTACTTTTGTTTTTGTTTTTTTTTTTTTTTTTTTACTATTATTAATAACAC')

    def test_molecule_pooling(self):
        for sample in ['AP1-P22-1-1_318','APKS2-P8-2-2_52']:
            for hd in [0,1]:

                f= pysam.AlignmentFile('./data/mini_nla_test.bam')
                it = singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    moleculeClass=singlecellmultiomics.molecule.Molecule,
                    fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
                    fragment_class_args={'umi_hamming_distance':hd}

                )
                molecule_count = 0
                for molecule in iter(it):
                    if molecule.get_sample()==sample:
                        molecule_count+=1
                if hd==0:
                    # it has one fragment with a sequencing error in the UMI
                    self.assertEqual(molecule_count,2)
                else:
                    self.assertEqual(molecule_count,1)



if __name__ == '__main__':
    unittest.main()
