#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment

"""
These tests check if the Molecule module is working correctly
"""

class TestMolecule(unittest.TestCase):

    def test_consensus(self):
        it = singlecellmultiomics.molecule.MoleculeIterator(
        alignments=f,
        moleculeClass=singlecellmultiomics.molecule.Molecule,
        fragmentClass=singlecellmultiomics.fragment.Fragment)

        for molecule in iter(it):
            #print(molecule.get_sample())
            if  molecule.get_sample()=='APKS3-P19-1-1_91':
                break

        assertEqual(''.join( list(molecule.get_consensus().values()) ), 'AGTTAGATATGGACTCTTCTTCAGACACTTTGTTTAAATTTTAAATTTTTTTCTGATTGCATATTACTAAAAATGTGTTATGAATATTTTCCATATCATTAAACATTCTTCTCAAGCATAACTTTAAATAACTATAGAAAATTTACGCTACTTTTGTTTTTGTTTTTTTTTTTTTTTTTTTACTATTATTAATAACAC')


if __name__ == '__main__':
    unittest.main()
