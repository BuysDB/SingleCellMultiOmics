#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from singlecellmultiomics.fastqProcessing.fastqIterator import FastqRecord
from singlecellmultiomics.barcodeFileParser.barcodeFileParser import BarcodeParser
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.CELSeq2 import CELSeq2_c8_u6_NH
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.scCHIC import SCCHIC_384w_c8_u3_cs2
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod
import pkg_resources

class TestUmiBarcodeDemux(unittest.TestCase):

    def test_UmiBarcodeDemuxMethod_matching_barcode(self):

        barcode_folder = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/barcodes/')
        index_folder = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/indices/')
        barcode_parser = BarcodeParser(barcode_folder)
        index_parser = BarcodeParser(index_folder)

        r1 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 1:N:0:GTGAAA',
          'ATCACACACTATAGTCATTCAGGAGCAGGTTCTTCAGGTTCCCTGTAGTTGTGTGGTTTTGAGTGAGTTTTTTAAT',
          '+',
          'AAAAA#EEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEE/EEEEEEEEEEEEEEEEEE'
        )
        r2 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 2:N:0:GTGAAA',
          'ACCCCAGATCAACGTTGGACNTCNNCNTTNTNCTCNGCACCNNNNCNNNCTTATNCNNNANNNNNNNNNNTNNGN',
          '+',
          '6AAAAEEAEE/AEEEEEEEE#EE##<#6E#A#EEE#EAEEA####A###EE6EE#E###E##########E##A#'
        )
        demux = UmiBarcodeDemuxMethod(umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            barcodeFileParser=barcode_parser,
            barcodeFileAlias='maya_384NLA',
            indexFileParser=index_parser,
            indexFileAlias='illumina_merged_ThruPlex48S_RP',
            random_primer_read=None,
            random_primer_length=6)

        demultiplexed_record = demux.demultiplex([r1,r2])
        # The barcode sequence is ACACACTA (first barcode)
        self.assertEqual( demultiplexed_record[0].tags['BC'], 'ACACACTA')
        self.assertEqual( demultiplexed_record[0].tags['bi'], 1) # 1 from version 0.1.12


    def test_CS2_NH_matching_barcode(self):

        barcode_folder = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/barcodes/')
        index_folder = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/indices/')
        barcode_parser = BarcodeParser(barcode_folder)
        index_parser = BarcodeParser(index_folder)

        seq = 'TATGAGCAATCACACACTATAGTCATTCAGGAGCAGGTTCTTCAGGTTCCCTGTAGTTGTGT'
        r1 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 1:N:0:GTGAAA',
          f'ATAATATCTGGGCA{seq}',
          '+',
          'AAAAA#EEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEE/EEEEEEEEEEEEEEEEEE'
        )
        r2 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 2:N:0:GTGAAA',
          'ACCCCAGATCAACGTTGGACNTCNNCNTTNTNCTCNGCACCNNNNCNNNCTTATNCNNNANNNNNNNNNNTNNGN',
          '+',
          '6AAAAEEAEE/AEEEEEEEE#EE##<#6E#A#EEE#EAEEA####A###EE6EE#E###E##########E##A#'
        )
        demux = CELSeq2_c8_u6_NH(
            barcodeFileParser=barcode_parser,
            indexFileParser=index_parser)

        demultiplexed_record = demux.demultiplex([r1,r2])
        # The barcode sequence is ACACACTA (first barcode)
        self.assertEqual( demultiplexed_record[0].tags['BC'], 'TCTGGGCA')
        self.assertEqual( demultiplexed_record[0].tags['bi'], 55)
        self.assertEqual( demultiplexed_record[0].tags['RX'], 'ATAATA')
        self.assertEqual( demultiplexed_record[0].sequence, seq)

    def test_CHICT(self):

        barcode_folder = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/barcodes/')
        index_folder = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/indices/')
        barcode_parser = BarcodeParser(barcode_folder)
        index_parser = BarcodeParser(index_folder)

        seq = 'TATGAGCAATCACACACTATAGTCATTCAGGAGCAGGTTCTTCAGGTTCCCTGTAGTTGTGT'
        r1 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 1:N:0:GTGAAA',
          f'ATAATATCTGGGCA{seq}',
          '+',
          'AAAAA#EEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEE/EEEEEEEEEEEEEEEEEE'
        )
        r2 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 2:N:0:GTGAAA',
          'ACCCCAGATCAACGTTGGACNTCNNCNTTNTNCTCNGCACCNNNNCNNNCTTATNCNNNANNNNNNNNNNTNNGN',
          '+',
          '6AAAAEEAEE/AEEEEEEEE#EE##<#6E#A#EEE#EAEEA####A###EE6EE#E###E##########E##A#'
        )
        demux = SCCHIC_384w_c8_u3_cs2(
            barcodeFileParser=barcode_parser,
            indexFileParser=index_parser)

        demultiplexed_record = demux.demultiplex([r1,r2])
        # The barcode sequence is ACACACTA (first barcode)
        self.assertEqual( demultiplexed_record[0].tags['BC'], 'TCTGGGCA')
        self.assertEqual( demultiplexed_record[0].tags['bi'], 55)
        self.assertEqual( demultiplexed_record[0].tags['MX'], 'CS2C8U6')
        self.assertEqual( demultiplexed_record[0].tags['RX'], 'ATAATA')
        self.assertEqual( demultiplexed_record[0].sequence, seq)



        seq = 'TATGAGCAATCACACACTATAGTCATTCAGGAGCAGGTTCTTCAGGTTCCCTGTAGTTGTGT'

        R1_seq = f'AAAAGAGCGCGT{seq}'
        r1 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 1:N:0:GTGAAA',
          R1_seq,
          '+',
          'E'*len(R1_seq)
        )
        r2 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 2:N:0:GTGAAA',
          'ACCCCAGATCAACGTTGGACNTCNNCNTTNTNCTCNGCACCNNNNCNNNCTTATNCNNNANNNNNNNNNNTNNGN',
          '+',
          '6AAAAEEAEE/AEEEEEEEE#EE##<#6E#A#EEE#EAEEA####A###EE6EE#E###E##########E##A#'
        )
        demux = SCCHIC_384w_c8_u3_cs2(
            barcodeFileParser=barcode_parser,
            indexFileParser=index_parser)

        demultiplexed_record = demux.demultiplex([r1,r2])
        # The barcode sequence is ACACACTA (first barcode)
        self.assertEqual( demultiplexed_record[0].tags['BC'], 'AGAGCGCG')
        self.assertEqual( demultiplexed_record[0].tags['bi'], 26)
        self.assertEqual( demultiplexed_record[0].tags['MX'], 'scCHIC384C8U3')
        self.assertEqual( demultiplexed_record[0].tags['RX'], 'AAA')
        self.assertEqual( demultiplexed_record[0].sequence, seq)


    def test_3DEC_UmiBarcodeDemuxMethod_matching_barcode(self):

        barcode_folder = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/barcodes/')
        barcode_parser = BarcodeParser(barcode_folder)

        r1 = FastqRecord(
          '@Cluster_s_1_1101_1000',
          'ATCACACACTATAGTCATTCAGGAGCAGGTTCTTCAGGTTCCCTGTAGTTGTGTGGTTTTGAGTGAGTTTTTTAAT',
          '+',
          'AAAAA#EEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEE/EEEEEEEEEEEEEEEEEE'
        )
        r2 = FastqRecord(
          '@Cluster_s_1_1101_1002',
          'ACCCCAGATCAACGTTGGACNTCNNCNTTNTNCTCNGCACCNNNNCNNNCTTATNCNNNANNNNNNNNNNTNNGN',
          '+',
          '6AAAAEEAEE/AEEEEEEEE#EE##<#6E#A#EEE#EAEEA####A###EE6EE#E###E##########E##A#'
        )
        demux = UmiBarcodeDemuxMethod(umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            barcodeFileParser=barcode_parser,
            barcodeFileAlias='maya_384NLA',
            indexFileParser=None,
            indexFileAlias='illumina_merged_ThruPlex48S_RP',
            random_primer_read=None,
            random_primer_length=6)

        demultiplexed_record = demux.demultiplex([r1,r2])
        # The barcode sequence is ACACACTA (first barcode)
        self.assertEqual( demultiplexed_record[0].tags['BC'], 'ACACACTA')
        self.assertEqual( demultiplexed_record[0].tags['bi'], 1)


if __name__ == '__main__':
    unittest.main()
