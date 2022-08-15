#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
import pysamiterators.iterators
import os

from singlecellmultiomics.molecule import MoleculeIterator, CHICMolecule
from singlecellmultiomics.fragment import CHICFragment

"""
These tests check if the Molecule module is working correctly
"""

class TestMolecule(unittest.TestCase):

    def test_chic_cigar_dedup(self):
        i = 0
        with pysam.AlignmentFile('./data/chic_test_region.bam') as alignments:

            for molecule in MoleculeIterator(alignments,CHICMolecule, CHICFragment):
                i+=1

        self.assertEqual(i,1)


    def test_eq(self):


        def get_chic_read(header,  qname, contig ='chr1', start=100, sequence='ATCGGG', cigar=None, umi='CAT', sample='CELL_1', is_reverse=False, read1=True, paired=False,proper_pair=True ):
                read = pysam.AlignedSegment(header)
                read.set_tag('SM',sample) # The sample to which the sample belongs is extracted from the SM tag
                read.set_tag('RX',umi) # The UMI is extracted from the RX tag
                read.set_tag('MX','scCHIC')
                # By default the molecule assignment is done based on the mapping location of read 1:
                read.reference_name = contig
                read.reference_start = start
                read.query_name = qname
                read.query_sequence = sequence
                read.is_reverse = is_reverse
                read.cigarstring = f'{len(sequence)}M' if cigar is None else cigar
                if read1:
                    read.is_read1 = True
                    read.is_read2 = False
                else:
                    read.is_read1 = False
                    read.is_read2 = True
                if paired:
                    read.is_paired = True
                    read.is_proper_pair = proper_pair
                return read

        from singlecellmultiomics.molecule import Molecule
        from singlecellmultiomics.fragment import Fragment
        import pysam
        # Create sam file to write some reads to:
        with pysam.AlignmentFile('test.sam','w',reference_names=['chr1','chr2'],reference_lengths=[1000,1000]) as test_sam:

            read_A = get_chic_read(test_sam.header, 'read_A')
            read_B = get_chic_read(test_sam.header, 'read_B', start=102)
            read_C = get_chic_read(test_sam.header, 'read_C', start=110)
            # A reverse read
            read_D = get_chic_read(test_sam.header, 'read_D', start=110-5, is_reverse=True)
            # A different sample read
            read_E = get_chic_read(test_sam.header, 'read_E', start=110-5, is_reverse=True, sample='CELL_2')

            # A umi error read
            read_F = get_chic_read(test_sam.header, 'read_F', start=110-5, is_reverse=True, sample='CELL_2', umi='CAG')

            # A proper pair
            read_G_a = get_chic_read(test_sam.header, 'read_G', start=102,  paired=True)
            read_G_b = get_chic_read(test_sam.header, 'read_G', start=130, read1=False, paired=True, is_reverse=True)

            # A stale R2
            read_H = get_chic_read(test_sam.header, 'read_H', start=130,  paired=True, proper_pair=False, read1=False)

            # Read softclipped from start
            read_I = get_chic_read(test_sam.header, 'soft_clipped_read', start=104, cigar='2S4M')

            frag_A = Fragment([read_A],umi_hamming_distance=0,assignment_radius=2)
            frag_B = Fragment([read_B],umi_hamming_distance=0,assignment_radius=2)
            frag_C = Fragment([read_C],umi_hamming_distance=0,assignment_radius=2)
            frag_D = Fragment([read_D],umi_hamming_distance=0,assignment_radius=2)
            frag_E = Fragment([read_E],umi_hamming_distance=0,assignment_radius=2)
            frag_F = Fragment([read_F],umi_hamming_distance=0,assignment_radius=2)
            frag_G = Fragment([read_G_a, read_G_b],umi_hamming_distance=0,assignment_radius=2)

            self.assertTrue(frag_A==frag_B)
            self.assertTrue(frag_A==frag_G)
            self.assertFalse(frag_A==frag_C)
            self.assertFalse(frag_A==frag_D)
            self.assertFalse(frag_C==frag_D)
            self.assertFalse(frag_E==frag_D)

            frag_A = CHICFragment([read_A],umi_hamming_distance=0,assignment_radius=2)
            frag_B = CHICFragment([read_B],umi_hamming_distance=0,assignment_radius=2)
            frag_C = CHICFragment([read_C],umi_hamming_distance=0,assignment_radius=2)
            frag_D = CHICFragment([read_D],umi_hamming_distance=0,assignment_radius=2)
            frag_E = CHICFragment([read_E],umi_hamming_distance=1,assignment_radius=2)
            frag_F = CHICFragment([read_F],umi_hamming_distance=1,assignment_radius=2)
            frag_G = CHICFragment([read_G_a, read_G_b],umi_hamming_distance=0,assignment_radius=2)
            frag_I = CHICFragment([read_I],umi_hamming_distance=1,assignment_radius=2)


            frag_H = CHICFragment([None, read_H],umi_hamming_distance=1,assignment_radius=2)
            self.assertFalse(frag_H.is_valid())
            self.assertTrue(frag_A.is_valid())
            self.assertTrue(frag_B.is_valid())

            self.assertTrue(frag_A==frag_B)

            self.assertEqual(frag_G, frag_I)
            self.assertFalse(frag_A==frag_C)
            self.assertFalse(frag_A==frag_D)
            self.assertFalse(frag_C==frag_D)
            self.assertFalse(frag_E==frag_D)


            self.assertTrue(frag_A==frag_B)
            self.assertFalse(frag_A==frag_C)

            molecule_A = CHICMolecule([frag_A])
            print(molecule_A.site_location[1])
            self.assertEqual(molecule_A.site_location[1],98)
            self.assertTrue(molecule_A==frag_B)
            self.assertFalse(molecule_A==frag_C)

            # Test construction of molecule:
            self.assertTrue( molecule_A.add_fragment(frag_B) )
            self.assertEqual(len(molecule_A),2)

            # Frag C cannot be added:
            self.assertFalse( molecule_A.add_fragment(frag_C) )
            self.assertEqual(len(molecule_A),2)

            # Frag D cannot be added:
            self.assertFalse( molecule_A.add_fragment(frag_D) )
            self.assertEqual(len(molecule_A),2)

            # Test moving of site location by read which is aligned more to left:
            molecule = CHICMolecule([frag_B])
            self.assertEqual(molecule.site_location[1],100)
            self.assertTrue( molecule.add_fragment(frag_A) )
            self.assertEqual(molecule.site_location[1],98)

            # Test moving of site location by read which is aligned more to right: (Reverse strand)
            molecule = CHICMolecule([frag_B])
            self.assertEqual(molecule.site_location[1],100)
            self.assertTrue( molecule.add_fragment(frag_A) )
            self.assertEqual(molecule.site_location[1],98)

            # Test umi distance:
            molecule = CHICMolecule([frag_E])
            self.assertTrue( molecule.add_fragment(frag_F) )
            read_F.set_tag('RX','GGG')
            frag_G = CHICFragment([read_F],umi_hamming_distance=1,assignment_radius=2)
            self.assertFalse( molecule.add_fragment(frag_G) )


            # Perform iteration
            reads = [read_A, read_B, read_C, read_D, read_E, read_F, (read_G_a, read_G_b), (None,read_H), (read_I,None)]
            molecules = list(MoleculeIterator( reads,
                            molecule_class = CHICMolecule,
                            fragment_class = CHICFragment,
                            yield_invalid=True,
                            fragment_class_args={
                                'assignment_radius':2,
                                'umi_hamming_distance':1
                            }))

            recovered_reads = []
            for molecule in molecules:
                for read in molecule.iter_reads():
                    recovered_reads.append(read)
                    print(read.query_name)
            self.assertEqual(len(recovered_reads),10)

        os.remove('test.sam')



    def test_chic_nocigar_dedup(self):
        i = 0
        with pysam.AlignmentFile('./data/chic_test_region.bam') as alignments:
            for molecule in MoleculeIterator(alignments,CHICMolecule, CHICFragment,fragment_class_args={'no_umi_cigar_processing':True}):
                i+=1
        self.assertEqual(i,2)

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
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.Fragment
            )
            pass

    def test_every_fragment_as_molecule(self):
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for i,m in enumerate(singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.Fragment,
            every_fragment_as_molecule=True
            )):
                pass
            self.assertEqual(i,340)

    def test_every_fragment_as_molecule_np_iterator(self):

        to_be_found = set()
        amount_of_r1s_to_be_found=0
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for read in f:
                if read.is_read1 and not read.is_read2:
                    to_be_found.add(read.query_name)
                    amount_of_r1s_to_be_found+=1
        found=set()
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            frag_count=0
            for i,m in enumerate(singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.Fragment,
            every_fragment_as_molecule=True,
            yield_invalid=True,
            iterator_class=pysamiterators.MatePairIteratorIncludingNonProper
            )):
                if m[0].has_R1():
                    frag_count+=1
                    found.add( m[0][0].query_name )


            diff = to_be_found.difference( found )
            tm = found.difference(to_be_found)
            print('missed:')
            print(diff)
            print('too much:')
            print(tm)
            self.assertEqual( len(diff), 0)
            self.assertEqual(frag_count,amount_of_r1s_to_be_found)
            self.assertEqual(frag_count,293)





    def test_Molecule_repr_stability(self):
        """Test if the molecule representation function is stable"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.Fragment
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

            # Test tags:
            a = molecules[0]
            a.write_tags()
            self.assertEqual(a[0][0].get_tag('TF') , 2)

        os.remove('test.sam')
    def test_NLA_Molecule_repr_stability(self):
        """Test if the molecule representation function is stable"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment
            )
            for molecule in it:
                str(molecule)

    def test_NLA_Molecule_ref_id(self):
        """Test if the molecule representation function is stable"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment
            )
            for molecule in it:
                self.assertEqual( molecule.get_a_reference_id(),0)



    def test_NLA_molecule_iterator_stability(self):
        """Test if the simplest molecule iterator works"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment
            )
            pass

    def test_consensus(self):
        """Test if a right consensus sequence can be produced from a noisy molecule"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
            fragment_class_args={
                'R1_primer_length':0,
                'R2_primer_length':6,
            }
            )
            for molecule in iter(it):
                #print(molecule.get_sample())
                if  molecule.get_sample()=='APKS3-P19-1-1_91':
                    break

            self.assertEqual(''.join( list(molecule.get_consensus().values()) ), 'CATGAGTTAGATATGGACTCTTCTTCAGACACTTTGTTTAAATTTTAAATTTTTTTCTGATTGCATATTACTAAAAATGTGTTATGAATATTTTCCATATCATTAAACATTCTTCTCAAGCATAACTTTAAATAACTGCATTATAGAAAATTTACGCTACTTTTGTTTTTGTTTTTTTTTTTTTTTTTTTACTATTATTAATAACACGGTGG')

    def test_fragment_sizes(self):

        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                alignments=f,
                molecule_class=singlecellmultiomics.molecule.Molecule,
                fragment_class_args={
                    'R1_primer_length':4,
                    'R2_primer_length':6,
                },

                fragment_class=singlecellmultiomics.fragment.NlaIIIFragment):

                # This is a dovetailed molecule, both R1 and R2 overshoot into the adapter
                if molecule.sample=='APKS2-P18-2-1_66' and molecule.umi=='CGC':
                    # The non-dovetailed region of this molecule is
                    # 164834866-164834975
                    # This means a fragment size of 109,
                    # the 4 bases of the CATG are not counted.
                    # the 6 bases of the random primer are also not counted
                    # Resulting in a fragment size of 109 - 10 = 101

                    start, end = molecule[0].get_safe_span()
                    self.assertEqual(
                        abs(end-start)
                        , 107)
                    self.assertEqual(
                        molecule.get_safely_aligned_length()
                        , 107)


    def test_get_match_mismatch_frequency(self):
        """Test if the matches and mismatches of reads in a molecule are counted properly"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class_args={
                'R1_primer_length':0,
                'R2_primer_length':0,
            },

            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment)
            for molecule in iter(it):
                #print(molecule.get_sample())
                if  molecule.get_sample()=='APKS3-P19-1-1_91':
                    break
        #print(molecule.get_match_mismatch_frequency())
        self.assertEqual( (628, 13), molecule.get_match_mismatch_frequency() )

    def test_get_consensus_gc_ratio(self):
        """Test if the gc ratio of a molecule is properly determined"""
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.Fragment,
            fragment_class_args={
                'R1_primer_length':0,
                'R2_primer_length':0,
            }
            )
            for molecule in iter(it):
                #print(molecule.get_sample())
                if  molecule.get_sample()=='APKS3-P19-1-1_91':
                    break

        self.assertAlmostEqual(0.23113207547169812, molecule.get_consensus_gc_ratio() )


    def test_deduplication(self):
        f= pysam.AlignmentFile('./data/mini_nla_test.bam')
        total_frag_count = 0
        for molecule in singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
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
                molecule_class=singlecellmultiomics.molecule.Molecule,
                fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
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

    def test_max_associated_fragments(self):

        for i in range(1,3):
            with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
                for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    molecule_class=singlecellmultiomics.molecule.Molecule,
                    fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
                    molecule_class_args={'max_associated_fragments':i},
                    fragment_class_args={'umi_hamming_distance':1}

                ):
                    self.assertTrue(len(molecule)<= i , f'{i}, {len(molecule)}')


    def test_rt_reaction_counting_HAMMING_1(self):
        # This is a dictionary containing molecules and the amount of rt reactions for location chr1:164834865

        f= pysam.AlignmentFile('./data/mini_nla_test.bam')
        it = singlecellmultiomics.molecule.MoleculeIterator(
            alignments=f,
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
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
            molecule_class=singlecellmultiomics.molecule.Molecule,
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
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
            molecule_class=singlecellmultiomics.molecule.AnnotatedNLAIIIMolecule,
            molecule_class_args={'features':features},
            fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
            fragment_class_args={'umi_hamming_distance':1}
        )

        hit_count = 0
        for molecule in it:
            molecule.annotate()
            if len(molecule.hits):
                hit_count+=1
        self.assertEqual(hit_count,2)

    """
    def test_classification_consensus(self):
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f, pysam.AlignmentFile('./data/consensus_write_test.bam','wb',header=f.header) as target_bam:
            molecule_iterator =  singlecellmultiomics.molecule.MoleculeIterator(
                alignments=f,
                molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
                fragment_class=singlecellmultiomics.fragment.NlaIIIFragment
            )

            classifier = singlecellmultiomics.molecule.train_consensus_model(
                        molecule_iterator,
                        mask_variants=None,
                        skip_already_covered_bases=False,
                        n_train=100)

            molecule_iterator =  singlecellmultiomics.molecule.MoleculeIterator(
                alignments=f,
                molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
                fragment_class=singlecellmultiomics.fragment.NlaIIIFragment
            )

            for i,molecule in enumerate(molecule_iterator):
                read_name=f'consensus_{i}'
                reads = molecule.deduplicate_to_single_CIGAR_spaced(target_bam,
                read_name, classifier,reference=None)

                read = molecule.deduplicate_to_single(target_bam)

                consensus = molecule.get_consensus(classifier)
                """
if __name__ == '__main__':
    unittest.main()
