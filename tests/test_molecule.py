#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
import pysamiterators.iterators
import os
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

    def test_molecule_iterator_stability(self):
        """Test if the simplest molecule iterator works"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.Fragment
            )
            pass

    def test_Molecule_repr_stability(self):
        """Test if the molecule representation function is stable"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.Fragment
            )
            for molecule in it:
                str(molecule)
    def test_iterable_molecule_iter(self):

        from singlecellmultiomics.molecule import MoleculeIterator
        from singlecellmultiomics.fragment import Fragment

        with pysam.AlignmentFile('test.sam','w',reference_names=['chr1','chr2'],reference_lengths=[1000,1000]) as test_sam:
            read_A = pysam.AlignedSegment(test_sam.header)
            read_A.set_tag('SM','CELL_1')
            read_A.set_tag('RX','CAT')
            read_A.reference_name = 'chr1'
            read_A.reference_start = 100
            read_A.query_sequence = 'ATCGGG'
            read_A.cigarstring = '6M'
            read_A.mapping_quality = 60

            read_B = pysam.AlignedSegment(test_sam.header)
            read_B.set_tag('SM','CELL_1')
            read_B.set_tag('RX','CAT')
            read_B.reference_name = 'chr1'
            read_B.reference_start = 100
            read_B.query_sequence = 'ATCGG'
            read_B.cigarstring = '5M'
            read_B.mapping_quality = 60

            read_C = pysam.AlignedSegment(test_sam.header)
            read_C.set_tag('SM','CELL_2')
            read_C.set_tag('RX','CAT')
            read_C.reference_name = 'chr1'
            read_C.reference_start = 100
            read_C.query_sequence = 'ATCGG'
            read_C.cigarstring = '5M'
            read_C.mapping_quality = 60

            reads = [  read_A,read_B,read_C ]
            mi = MoleculeIterator( reads , yield_invalid=True)
            molecules=[]
            for molecule in mi:
                molecules.append(molecule)

            self.assertEqual(len(molecules),2)
            self.assertEqual(max( (len(m) for m in molecules) ),2)
            self.assertEqual(min( (len(m) for m in molecules) ),1)

        os.remove('test.sam')
    def test_NLA_Molecule_repr_stability(self):
        """Test if the molecule representation function is stable"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment
            )
            for molecule in it:
                str(molecule)

    def test_NLA_Molecule_ref_id(self):
        """Test if the molecule representation function is stable"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment
            )
            for molecule in it:
                self.assertEqual( molecule.get_a_reference_id(),0)



    def test_NLA_molecule_iterator_stability(self):
        """Test if the simplest molecule iterator works"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment
            )
            pass

    def test_consensus(self):
        """Test if a right consensus sequence can be produced from a noisy molecule"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.Fragment,
            fragment_class_args={
                'R1_primer_length':4,
                'R2_primer_length':6,
            }
            )
            for molecule in iter(it):
                #print(molecule.get_sample())
                if  molecule.get_sample()=='APKS3-P19-1-1_91':
                    break

            self.assertEqual(''.join( list(molecule.get_consensus().values()) ), 'AGTTAGATATGGACTCTTCTTCAGACACTTTGTTTAAATTTTAAATTTTTTTCTGATTGCATATTACTAAAAATGTGTTATGAATATTTTCCATATCATTAAACATTCTTCTCAAGCATAACTTTAAATAACTATAGAAAATTTACGCTACTTTTGTTTTTGTTTTTTTTTTTTTTTTTTTACTATTATTAATAACAC')

    def test_fragment_sizes(self):

        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                alignments=f,
                moleculeClass=singlecellmultiomics.molecule.Molecule,
                fragment_class_args={
                    'R1_primer_length':4,
                    'R2_primer_length':6,
                },

                fragmentClass=singlecellmultiomics.fragment.Fragment):

                # This is a dovetailed molecule, both R1 and R2 overshoot into the adapter
                if molecule.sample=='APKS2-P18-2-1_66' and molecule.umi=='CGC':
                    # The non-dovetailed region of this molecule is
                    # 164834866-164834975
                    # This means a fragment size of 109,
                    # the 4 bases of the CATG are not counted.
                    # the 6 bases of the random primer are also not counted
                    # Resulting in a fragment size of 109 - 10 = 101
                    self.assertEqual(
                        abs(molecule[0].span[1]-molecule[0].span[2])
                        , 101)
                    self.assertEqual(
                        molecule.get_safely_aligned_length()
                        , 101)


    def test_get_match_mismatch_frequency(self):
        """Test if the matches and mismatches of reads in a molecule are counted properly"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragment_class_args={
                'R1_primer_length':4,
                'R2_primer_length':6,
            },

            fragmentClass=singlecellmultiomics.fragment.Fragment)
            for molecule in iter(it):
                #print(molecule.get_sample())
                if  molecule.get_sample()=='APKS3-P19-1-1_91':
                    break
        #print(molecule.get_match_mismatch_frequency())
        self.assertEqual( (583, 13), molecule.get_match_mismatch_frequency() )

    def test_get_consensus_gc_ratio(self):
        """Test if the gc ratio of a molecule is properly determined"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.Fragment,
            fragment_class_args={
                'R1_primer_length':4,
                'R2_primer_length':6,
            }
            )
            for molecule in iter(it):
                #print(molecule.get_sample())
                if  molecule.get_sample()=='APKS3-P19-1-1_91':
                    break

        self.assertAlmostEqual(0.2070707, molecule.get_consensus_gc_ratio() )


    def test_deduplication(self):
        f= pysam.AlignmentFile('./data/mini_nla_test.bam')
        total_frag_count = 0
        for molecule in singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
            fragment_class_args={'umi_hamming_distance':0},
            pooling_method=0,
            yield_invalid=True
        ):
            if 'APKS2-P14-1-1' in molecule.sample:
                for fragment in molecule:

                    is_dup = False
                    for read in fragment.reads:
                        if read is not None and read.is_duplicate:
                            is_dup = True
                    if not is_dup:
                        total_frag_count+=1
        self.assertEqual(total_frag_count,13)

    def _pool_test(self, pooling_method=0,hd=0):
        for sample in ['AP1-P22-1-1_318','APKS2-P8-2-2_52']:

            f= pysam.AlignmentFile('./data/mini_nla_test.bam')
            it = singlecellmultiomics.molecule.MoleculeIterator(
                alignments=f,
                moleculeClass=singlecellmultiomics.molecule.Molecule,
                fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
                fragment_class_args={'umi_hamming_distance':hd},
                pooling_method=pooling_method
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

    def test_molecule_pooling_vanilla_exact_umi(self):
        self._pool_test(0,0)

    def test_molecule_pooling_nlaIIIoptim_exact_umi(self):
        self._pool_test(1,0)

    def test_molecule_pooling_vanilla_umi_mismatch(self):
        self._pool_test(0,1)

    def test_molecule_pooling_nlaIIIoptim_umi_mismatch(self):
        self._pool_test(1,1)

    def test_rt_reaction_counting_HAMMING_1(self):
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
                obtained_rt_count[molecule.get_sample()]  = len( molecule.get_rt_reaction_fragment_sizes(max_N_distance=1) )

        # Validate:
        for sample, truth in hand_curated_truth.items():
            self.assertEqual( obtained_rt_count.get(sample,0),truth.get('rt_count',-1) )


    def test_rt_reaction_counting_HAMMING_0(self):
        # This is a dictionary containing molecules and the amount of rt reactions for location chr1:164834865

        f= pysam.AlignmentFile('./data/mini_nla_test.bam')
        it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
            fragment_class_args={'umi_hamming_distance':1}

        )

        hand_curated_truth = {
            # hd: 0:
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
                obtained_rt_count[molecule.get_sample()]  = len( molecule.get_rt_reaction_fragment_sizes(max_N_distance=0) )

        # Validate:
        for sample, truth in hand_curated_truth.items():
            self.assertEqual( obtained_rt_count.get(sample,0),truth.get('rt_count',-1) )


    def test_feature_molecule(self):
        import singlecellmultiomics.features
        features = singlecellmultiomics.features.FeatureContainer()
        features.addFeature(chromosome='chr1',start=164835366, end=164835535, data=('test4'),name='test4')
        features.sort()

        f= pysam.AlignmentFile('./data/mini_nla_test.bam')
        it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            moleculeClass=singlecellmultiomics.molecule.AnnotatedNLAIIIMolecule,
            molecule_class_args={'features':features},
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
            fragment_class_args={'umi_hamming_distance':1}
        )

        hit_count = 0
        for molecule in it:
            molecule.annotate()
            if len(molecule.hits):
                hit_count+=1
        self.assertEqual(hit_count,2)

    def test_classification_consensus(self):
        """This only tests if no error is raised"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f, pysam.AlignmentFile('./data/consensus_write_test.bam','wb',header=f.header) as target_bam:
            molecule_iterator =  singlecellmultiomics.molecule.MoleculeIterator(
                alignments=f,
                moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
                fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment
            )

            classifier = singlecellmultiomics.molecule.train_consensus_model(
                        molecule_iterator,
                        mask_variants=None,
                        n_train=100)

            molecule_iterator =  singlecellmultiomics.molecule.MoleculeIterator(
                alignments=f,
                moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
                fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment
            )

            for i,molecule in enumerate(molecule_iterator):
                read_name=f'consensus_{i}'
                reads = molecule.deduplicate_to_single_CIGAR_spaced(target_bam,
                read_name, classifier,reference=None)

                read = molecule.deduplicate_to_single(target_bam)

                consensus = molecule.get_consensus(classifier)

if __name__ == '__main__':
    unittest.main()
