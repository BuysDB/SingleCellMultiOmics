#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
import pysamiterators.iterators
import os
from singlecellmultiomics.fragment import Fragment, CHICFragment
from singlecellmultiomics.utils import create_MD_tag
from singlecellmultiomics.utils import complement
"""
These tests check if the Molecule module is working correctly
"""
class TestTAPs(unittest.TestCase):


    def test_random_priming(self):
        temp_folder = 'data'

        enable_ref_write=False
        ref_path = f'{temp_folder}/chic_ref.fa'
        alignments_path = f'{temp_folder}/chic_test_alignments.bam'

        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
        # Create reference bam file

        refseq = 'TTAATCATGAAACCGTGGAGGCAAATCGGAGTGTAAGGCTTGACTGGATTCCTACGTTGCGTAGGTTCATGGGGGG'
        if enable_ref_write:
            with open(ref_path, 'w') as f:
                f.write(f">chr1\n{refseq}\n>chr2\n{complement(refseq)}\n""")

            # This command needs to finish, which is not working properly during testing
            pysam.faidx(ref_path)

        # CATG at base 5
        # Create BAM file with NLA fragment:

        def get_reads():
            with pysam.AlignmentFile(alignments_path_unsorted,'wb',reference_names=['chr1'],reference_lengths=[len(refseq)]) as bam:


                ### Nla III mate pair example, containing 2 CpGs and 1 call on the wrong strand
                read_A = pysam.AlignedSegment(bam.header)
                read_A.reference_name = 'chr1'
                read_A.reference_start = 5
                # Before last A is a bogus G>A conversion to test strandness:
                read_A.query_sequence = 'CATGAAACCGTGGAGGCAAATTGGAGTAT'
                read_A.cigarstring = f'{len(read_A.query_sequence)}M'
                read_A.qual = 'A'*len(read_A.query_sequence)
                read_A.mapping_quality = 60
                read_A.query_name = 'EX1_GA_CONV_2x_CpG_TAPS'
                read_A.set_tag('SM', 'Cell_A')
                read_A.is_read1 = True
                read_A.is_read2 = False
                read_A.set_tag('lh','TG')
                # Set substitution tag:
                read_A.set_tag('MD',
                               create_MD_tag(
                                       refseq[read_A.reference_start:read_A.reference_end], read_A.query_sequence))
                read_A.is_paired = True
                read_A.is_proper_pair = True

                # Create a second read which is a mate of the previous
                read_B = pysam.AlignedSegment(bam.header)
                read_B.reference_name = 'chr1'
                read_B.reference_start = 25
                read_B.query_sequence = refseq[25:60].replace('TGT','TAT').replace('CG', 'TG')
                read_B.cigarstring = f'{len(read_B.query_sequence)}M'
                read_B.qual = 'A'*len(read_B.query_sequence)
                read_B.mapping_quality = 60
                read_B.is_read2 = True
                read_B.is_read1 = False
                read_B.is_reverse = True
                read_B.query_name = 'EX1_GA_CONV_2x_CpG_TAPS'
                read_B.set_tag('SM', 'Cell_A')
                read_B.set_tag('lh','TG')
                read_B.set_tag('MD',
                           create_MD_tag(refseq[read_B.reference_start:read_B.reference_end],
                                         read_B.query_sequence,
                                   ))
                read_B.is_paired = True
                read_B.is_proper_pair = True

                read_A.next_reference_id = read_B.reference_id
                read_A.next_reference_start = read_B.reference_start
                read_B.next_reference_id = read_A.reference_id
                read_B.next_reference_start = read_A.reference_start

                read_A.mate_is_reverse = read_B.is_reverse
                read_B.mate_is_reverse = read_A.is_reverse

                return read_A, read_B



        alignments_path_unsorted = f'{alignments_path}.unsorted.bam'

        read_A, read_B = get_reads()
        read_B.set_tag('rS','AATTAA') # Set random primer sequence
        ## Act on reads with random primer set:
        frag = Fragment([read_A, read_B])
        self.assertTrue(frag.R2_primer_length==0)
        self.assertTrue(frag.unsafe_trimmed)

        ## Act on reads without random primer
        read_A, read_B = get_reads()
        frag = Fragment([read_A, read_B], R2_primer_length=10)
        self.assertEqual(frag.R2_primer_length, 10)
        self.assertFalse(frag.unsafe_trimmed)

        read_A, read_B = get_reads()
        frag = CHICFragment([read_A, read_B], R2_primer_length=0)
        self.assertEqual(frag.R2_primer_length, 0)
        self.assertFalse(frag.unsafe_trimmed)

        read_A, read_B = get_reads()
        read_B.set_tag('rS','AATTAA') # Set random primer sequence
        frag = CHICFragment([read_A, read_B], R2_primer_length=0)
        self.assertEqual(frag.R2_primer_length, 0)
        self.assertTrue(frag.unsafe_trimmed)

        read_A, read_B = get_reads()
        read_B.set_tag('MX','scCHIC384C8U3l')
        frag = CHICFragment([read_A, read_B])
        self.assertEqual(frag.R2_primer_length, 0)
        self.assertFalse(frag.unsafe_trimmed)

        read_A, read_B = get_reads()
        read_B.set_tag('MX','scCHIC384C8U3')
        frag = CHICFragment([read_A, read_B])
        self.assertEqual(frag.R2_primer_length,6)
        self.assertFalse(frag.unsafe_trimmed)

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
