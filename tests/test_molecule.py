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

class TestMolecule(unittest.TestCase):

    def test_a_pysam_iterators(self):
        """Test if the pysamiterators package yields the proper amount or mate pairs"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for i,(R1,R2) in enumerate(pysamiterators.iterators.MatePairIterator(f)):
                pass
            self.assertEqual(i,224)

    def test_consensus(self):
        """Test if a right consensus sequence can be produced from a noisy molecule"""
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

    def test_rt_reaction_counting(self):
        # This is a dictionary containing molecules and the amount of rt reactions for location chr1:164834865

        f= pysam.AlignmentFile('./data/mini_nla_test.bam')
        it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
            fragment_class_args={'umi_hamming_distance':1}

        )

        hand_curated_truth = {
            # hd: 1:
            # hard example with N in random primer and N in UMI:
            'APKS2-P8-2-2_52':{'rt_count':3},

            # Simple examples:
            'APKS2-P18-1-1_318':{'rt_count':1},
            'APKS2-P18-1-1_369':{'rt_count':2},
            'APKS2-P18-2-1_66':{'rt_count':1},
            'APKS2-P18-2-1_76':{'rt_count':1},
            'APKS2-P18-2-1_76':{'rt_count':1},
        }


        obtained_rt_count = {}
        for molecule in it:
            site = molecule.get_cut_site()
            if site is not None and site[1]==164834865:
                obtained_rt_count[molecule.get_sample()]  = len( molecule.get_rt_reaction_fragment_sizes() )

        # Validate:
        for sample, truth in hand_curated_truth.items():
            self.assertEqual( obtained_rt_count.get(sample,0),truth.get('rt_count',-1) )


        #def test_rt_reaction_sizes(self):



if __name__ == '__main__':
    unittest.main()
