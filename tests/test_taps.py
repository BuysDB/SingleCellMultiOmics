#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import pysam
import os
import random

from singlecellmultiomics.molecule import TAPSNlaIIIMolecule, TAPS, TAPSCHICMolecule
from singlecellmultiomics.fragment import NlaIIIFragment,CHICFragment
from singlecellmultiomics.utils import create_MD_tag

class TestTAPs(unittest.TestCase):

    def test_all(self):
        # TAPS test cases

        ## Set up test dataset:
        temp_folder = 'tmp/scmo'
        ref_path = f'{temp_folder}/ref.fa'
        alignments_path = f'{temp_folder}/alignments.bam'
        tagged_path = f'{temp_folder}/alignments_tagged.bam'

        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
        # Create reference bam file
        refseq = 'TTAATCATGAAACCGTGGAGGCAAATCGGAGTGTAAGGCTTGACTGGATTCCTACGTTGCGTAGGTT'
        with open(ref_path, 'w') as f:
            f.write(f""">chr1
        {refseq}
        """)

        pysam.faidx(ref_path)

        # CATG at base 5
        # Create BAM file with NLA fragment:

        alignments_path_unsorted = f'{alignments_path}.unsorted.bam'
        with pysam.AlignmentFile(alignments_path_unsorted, 'wb', reference_names=['chr1'],
                                 reference_lengths=[100]) as bam:
            ### Nla III mate pair example, containing 2 CpGs and 1 call on the wrong strand
            read_A = pysam.AlignedSegment(bam.header)
            read_A.reference_name = 'chr1'
            read_A.reference_start = 5
            # Before last A is a bogus G>A conversion to test strandness:
            read_A.query_sequence = 'CATGAAACCGTGGAGGCAAATTGGAGTAT'
            read_A.cigarstring = f'{len(read_A.query_sequence)}M'
            read_A.qual = 'A' * len(read_A.query_sequence)
            read_A.mapping_quality = 60
            read_A.query_name = 'EX1_GA_CONV_2x_CpG_TAPS'
            read_A.set_tag('SM', 'Cell_A')
            read_A.is_read1 = True
            read_A.is_read2 = False
            read_A.set_tag('lh', 'TG')
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
            read_B.query_sequence = refseq[25:60].replace('TGT', 'TAT').replace('CG', 'TG')
            read_B.cigarstring = f'{len(read_B.query_sequence)}M'
            read_B.qual = 'A' * len(read_B.query_sequence)
            read_B.mapping_quality = 60
            read_B.is_read2 = True
            read_B.is_read1 = False
            read_B.is_reverse = True
            read_B.query_name = 'EX1_GA_CONV_2x_CpG_TAPS'
            read_B.set_tag('SM', 'Cell_A')
            read_B.set_tag('lh', 'TG')
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

            bam.write(read_A)
            bam.write(read_B)

            ### Nla III mate pair example, dove tailed over random primer
            # , containing 1 CpGs and one 1 call in the dove tail which should not be called
            read_C = pysam.AlignedSegment(bam.header)
            read_C.reference_name = 'chr1'
            read_C.reference_start = 5
            read_C.query_sequence = 'CATGAAACCGTGGAGGC'.replace('ACC', 'ATC').replace('AGGC', 'CGGT')
            read_C.cigarstring = f'{len(read_C.query_sequence)}M'
            read_C.qual = 'A' * len(read_C.query_sequence)
            read_C.mapping_quality = 60
            read_C.query_name = 'EX2_GA_DOVE'
            read_C.set_tag('SM', 'Cell_A')
            read_C.is_read1 = True
            read_C.set_tag('lh', 'TG')
            # Set substitution tag:
            read_C.set_tag('MD',
                           create_MD_tag(
                               refseq[read_C.reference_start:read_C.reference_end],
                               read_C.query_sequence))

            read_C.is_paired = True
            read_C.is_proper_pair = True

            # Create a second read which is a mate of the previous
            read_D = pysam.AlignedSegment(bam.header)
            read_D.reference_name = 'chr1'
            read_D.reference_start = 10
            read_D.query_sequence = refseq[10:15].replace('ACC', 'GTC')
            read_D.cigarstring = f'{len(read_D.query_sequence)}M'
            read_D.qual = 'A' * len(read_D.query_sequence)
            read_D.mapping_quality = 60
            read_D.is_read2 = True
            read_D.is_read1 = False
            read_D.is_reverse = True
            read_D.query_name = 'EX2_GA_DOVE'
            read_D.set_tag('SM', 'Cell_A')
            read_D.set_tag('lh', 'TG')
            read_D.set_tag('MD',
                           create_MD_tag(refseq[read_D.reference_start:read_D.reference_end],
                                         read_D.query_sequence,
                                         ))
            read_D.is_paired = True
            read_D.is_proper_pair = True

            read_C.next_reference_id = read_D.reference_id
            read_C.next_reference_start = read_D.reference_start
            read_D.next_reference_id = read_C.reference_id
            read_D.next_reference_start = read_C.reference_start

            bam.write(read_C)
            bam.write(read_D)

        pysam.sort(alignments_path_unsorted, '-o', alignments_path)
        pysam.index(alignments_path)

        ####

        taps = TAPS()
        reference = pysam.FastaFile(ref_path)

        molecule = TAPSNlaIIIMolecule(
            NlaIIIFragment([read_A, read_B], invert_strand=False),
            reference=reference,
            taps=taps,
            taps_strand='F'
        )
        molecule.__finalise__()

        calls = molecule.methylation_call_dict

        self.assertEqual( calls['chr1', 54]['context'], 'Z')
        self.assertEqual( calls['chr1', 26]['context'], 'Z')
        self.assertNotIn(  ('chr1', 26 + 6), calls)

        molecule = TAPSNlaIIIMolecule(
            NlaIIIFragment([read_C, read_D]),
            reference=reference,
            taps=taps,
            taps_strand='F'
        )
        molecule.__finalise__()

        calls = molecule.methylation_call_dict

        self.assertEqual(calls['chr1', 12]['context'], 'X')

        # Check that dove tail is not included:
        self.assertNotIn(('chr1', 21), calls)


if __name__ == '__main__':
    unittest.main()
