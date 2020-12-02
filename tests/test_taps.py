#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import pysam
import os

from singlecellmultiomics.molecule import TAPSNlaIIIMolecule, TAPS
from singlecellmultiomics.fragment import NlaIIIFragment
from singlecellmultiomics.utils import create_MD_tag

from singlecellmultiomics.utils import complement


class TestTAPs(unittest.TestCase):

    def test_all(self):
        temp_folder = 'data'

        enable_ref_write=True


        ref_path = f'{temp_folder}/ref.fa'
        alignments_path = f'{temp_folder}/alignments.bam'

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

        alignments_path_unsorted = f'{alignments_path}.unsorted.bam'
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

            bam.write(read_A)
            bam.write(read_B)

            ### Nla III mate pair example, dove tailed over random primer
            # , containing 1 CpGs and one 1 call in the dove tail which should not be called
            read_C = pysam.AlignedSegment(bam.header)
            read_C.reference_name = 'chr1'
            read_C.reference_start = 5
            read_C.query_sequence = 'CATGAAACCGTGGAGGC'.replace('ACC','ATC').replace('AGGC','CGGT')
            read_C.cigarstring = f'{len(read_C.query_sequence)}M'
            read_C.qual = 'A'*len(read_C.query_sequence)
            read_C.mapping_quality = 60
            read_C.query_name = 'EX2_GA_DOVE'
            read_C.set_tag('SM', 'Cell_A')
            read_C.is_read1 = True
            read_C.set_tag('lh','TG')
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
            read_D.query_sequence = refseq[10:15].replace('ACC','GTC')
            read_D.cigarstring = f'{len(read_D.query_sequence)}M'
            read_D.qual = 'A'*len(read_D.query_sequence)
            read_D.mapping_quality = 60
            read_D.is_read2 = True
            read_D.is_read1 = False
            read_D.is_reverse = True
            read_D.query_name = 'EX2_GA_DOVE'
            read_D.set_tag('SM', 'Cell_A')
            read_D.set_tag('lh','TG')
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

            read_C.mate_is_reverse = read_D.is_reverse
            read_D.mate_is_reverse = read_C.is_reverse

            bam.write(read_C)
            bam.write(read_D)

            ########################################
            # Reverse dovetailed (2 way) alignment #
            ########################################

            read_E = pysam.AlignedSegment(bam.header)
            read_E.reference_name = 'chr1'
            read_E.query_sequence = refseq[2:71].replace('CATGAA','CATAAA').replace('CGG','CAG')
            read_E.reference_start = 71 - len(read_E.query_sequence)
            read_E.cigarstring = f'{len(read_E.query_sequence)}M'
            read_E.qual = 'A'*len(read_E.query_sequence)
            read_E.mapping_quality = 60
            read_E.query_name = 'EX2_GA_2xDOVE_rev'
            read_E.set_tag('SM', 'Cell_A')
            read_E.is_read2 = False
            read_E.is_read1 = True
            read_E.set_tag('lh','TG')
            read_E.is_reverse = True
            # Set substitution tag:
            read_E.set_tag('MD',
                           create_MD_tag(
                                   refseq[read_E.reference_start:read_E.reference_end],
                               read_E.query_sequence))
            read_E.set_tag('ri','read_E')
            read_E.is_paired = True
            read_E.is_proper_pair = True


            # Create a second read which is a mate of the previous
            read_F = pysam.AlignedSegment(bam.header)
            read_F.reference_name = 'chr1'
            read_F.reference_start = 10
            read_F.query_sequence = refseq[10:74].replace('CGG','CAG').replace('GGGG','GAGG')
            read_F.cigarstring = f'{len(read_F.query_sequence)}M'
            read_F.qual = 'A'*len(read_F.query_sequence)
            read_F.mapping_quality = 60
            read_F.is_read1 = False
            read_F.is_read2 = True
            read_F.is_reverse = False
            read_F.query_name = 'EX2_GA_2xDOVE_rev'
            read_F.set_tag('ri','read_F')
            read_F.set_tag('SM', 'Cell_A')
            read_F.set_tag('lh','TG')
            read_F.set_tag('MD',
                       create_MD_tag(refseq[read_F.reference_start:read_F.reference_end],
                                     read_F.query_sequence,
                               ))
            read_F.is_paired = True
            read_F.is_proper_pair = True


            read_F.mate_is_reverse = read_E.is_reverse
            read_E.mate_is_reverse = read_F.is_reverse

            read_E.next_reference_id = read_F.reference_id
            read_E.next_reference_start = read_F.reference_start
            read_F.next_reference_id = read_E.reference_id
            read_F.next_reference_start = read_E.reference_start

            bam.write(read_E)
            bam.write(read_F)


        pysam.sort(alignments_path_unsorted, '-o', alignments_path)
        pysam.index(alignments_path)

        taps = TAPS()
        with pysam.FastaFile(ref_path) as reference:

            self.assertEqual(reference.fetch('chr1', 26, 26 + 3),'CGG')
            molecule = TAPSNlaIIIMolecule(
                NlaIIIFragment([read_A, read_B]),
                reference=reference,
                taps=taps,
                taps_strand='F'
            )
            molecule.__finalise__()

            calls = molecule.methylation_call_dict
            print(calls)
            print(calls[('chr1', 54)])
            self.assertEqual( calls['chr1', 54]['context'], 'Z')
            self.assertEqual( calls['chr1', 26]['context'], 'Z')
            self.assertNotIn(  ('chr1', 26 + 6), calls)


            molecule = TAPSNlaIIIMolecule(
                NlaIIIFragment([read_E, read_F]),
                reference =reference,
                taps = taps,
                taps_strand='F'
            )
            molecule.__finalise__()

            # Test dove-tail detection:
            self.assertNotIn( ('chr1', 71) , molecule.methylation_call_dict)
            self.assertNotIn(('chr1', 8) , molecule.methylation_call_dict)



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
