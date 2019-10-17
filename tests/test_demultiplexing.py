#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from singlecellmultiomics.fastqProcessing.fastqIterator import FastqRecord
from singlecellmultiomics.barcodeFileParser.barcodeFileParser import BarcodeParser
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
        self.assertEqual( demultiplexed_record[0].tags['BI'], 0)


if __name__ == '__main__':
    unittest.main()
