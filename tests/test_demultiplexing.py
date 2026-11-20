#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from singlecellmultiomics.fastqProcessing.fastqIterator import FastqRecord
from singlecellmultiomics.barcodeFileParser.barcodeFileParser import BarcodeParser
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.CELSeq2 import CELSeq2_c8_u6_NH
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.scCHIC import SCCHIC_384w_c8_u3_cs2
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.DamID import DamID2andT_SCA,DamID2_SCA
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import TaggedRecord, UmiBarcodeDemuxMethod, TagDefinitions
from singlecellmultiomics.modularDemultiplexer.demultiplexingStrategyLoader import DemultiplexingStrategyLoader
from singlecellmultiomics.utils import reverse_complement
import importlib.resources
import subprocess
import sys
import tempfile
import os
class Dummy_pysam_read:
    def __init__(self, header:str):
        self.query_name = header

class TestUmiBarcodeDemux(unittest.TestCase):

    def test_UmiBarcodeDemuxMethod_matching_barcode(self):

        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        index_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/'))
        barcode_parser = BarcodeParser(barcode_folder, lazyLoad='*')
        index_parser = BarcodeParser(index_folder, lazyLoad='*')

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
        self.assertEqual( demultiplexed_record[0].tags['Is'], 'NS500414' ) # Instrument should just match NS500414 (no preceding @)



    def test_CS2_NH_matching_barcode(self):

        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        index_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/'))
        barcode_parser = BarcodeParser(barcode_folder, lazyLoad='*')
        index_parser = BarcodeParser(index_folder, lazyLoad='*')

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

    def construct_tchic_read(self,crx,ccb,trx,tcb,mr,linker):
        seq = f'{crx}{ccb}{linker}{trx}{tcb}TTTTTTTTTTTTTTTTTTTTT{mr}'
        r1 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 1:N:0:GTGAAA',
          f'{seq}',
          '+',
          'A'*len(seq)
        )
        r2 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 2:N:0:GTGAAA',
          reverse_complement(seq),
          '+',
          'A'*len(seq)
        )
        return r1, r2

    def test_TCHIC(self):

        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        index_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/'))
        barcode_parser = BarcodeParser(barcode_folder, lazyLoad='*',)
        index_parser = BarcodeParser(index_folder, lazyLoad='*')

        crx = 'TAT'
        ccb = 'TAAGTGCT'
        trx = 'CTGTTG'
        tcb = 'ACAGAAGC'
        mr = 'TGAGAGAGAGAGAGAGAGAGAGAGC'
        linker = 'TATC'
        r1,r2 = self.construct_tchic_read(crx,ccb,trx,tcb,mr,linker)

        demux = SCCHIC_384w_c8_u3_cs2(
            barcodeFileParser=barcode_parser,
            indexFileParser=index_parser)

        demultiplexed_record = demux.demultiplex([r1,r2])
        # The barcode sequence is ACACACTA (first barcode)
        self.assertEqual( demultiplexed_record[0].tags['BC'], ccb)
        self.assertEqual( demultiplexed_record[0].tags['bi'], 225)
        self.assertEqual( demultiplexed_record[0].tags['dt'], 'VASA')
        self.assertEqual( demultiplexed_record[0].tags['RX'], crx)
        self.assertEqual( demultiplexed_record[0].tags['rx'], trx)
        self.assertEqual( demultiplexed_record[1].sequence, reverse_complement(mr)[:len(mr)-4])


        r1 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 1:N:0:GTGAAA',
          'GGCGACGTCCTTCACTATAGGGAGTTCTACAGTTCGACGATCCTTAAATGGTGAGTTTTTTTTTTTTTTTTTTTTTTTGACCGACGGTCCCCCCGGGACCC',
          '+',
          'A'*len('GGCGACGTCCTTCACTATAGGGAGTTCTACAGTTCGACGATCCTTAAATGGTGAGTTTTTTTTTTTTTTTTTTTTTTTGACCGACGGTCCCCCCGGGACCC')
        )
        r2 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 2:N:0:GTGAAA',
          'CGATCCTTAAATGGTGAGTTTTTTTTTTTTTTTTTTTTTTTGACCGACGGTCCCCCCGGGACCCGACGGCGCGACGACGCCCGGGGCGCACTGGGGACAGT',
          '+',
          'A'*len('CGATCCTTAAATGGTGAGTTTTTTTTTTTTTTTTTTTTTTTGACCGACGGTCCCCCCGGGACCCGACGGCGCGACGACGCCCGGGGCGCACTGGGGACAGT')
        )
        demultiplexed_record = demux.demultiplex([r1,r2])
        self.assertEqual( demultiplexed_record[0].tags['BC'], 'GACGTCCT')
        self.assertEqual( demultiplexed_record[0].tags['bi'], 214)
        self.assertEqual( demultiplexed_record[0].tags['dt'], 'VASA')
        self.assertEqual( demultiplexed_record[0].tags['RX'], 'GGC')
        self.assertEqual( demultiplexed_record[0].tags['rx'], 'CTTAAA')
        assert 'oh' not in demultiplexed_record[0].tags # The presence of the oh tag means the record was not parsed as Illumina data


    def construct_read_pair(self, prefix, content):
        seq = f'{prefix}{content}'
        r1 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 1:N:0:GTGAAA',
          f'{seq}',
          '+',
          'A'*len(seq)
        )
        r2 = FastqRecord(
          '@NS500414:628:H7YVNBGXC:1:11101:15963:1046 2:N:0:GTGAAA',
          reverse_complement(seq),
          '+',
          'A'*len(seq)
        )
        return r1, r2
    
    def test_decode_tagged(self):
        tr = TaggedRecord(TagDefinitions)

        tr.fromTaggedBamRecord(Dummy_pysam_read("Is:@LH00371;RN:377;Fc:23JCYHLT3;La:1;Ti:2146;CX:50248;CY:15716;Fi:N;CN:0;aa:GACGAC;oh:LH00371:377:23JCYHLT3:1:2146:50248:15716;LY:JvB-203-RPE1-scFP-EdU-seq-super-pl17;RX:AGC;RQ:jyy;bi:126;bc:GCGCTACG;MX:scCHIC384C8U3;BC:GCGCTACG;rS:CCTAGC;lh:TA;lq:OO"))
        assert tr.tags['RN'] == '377'
        assert tr.tags['Is'] == 'LH00371'
        assert tr.tags['RX'] == 'AGC'
        assert tr.tags['oh'] == 'LH00371:377:23JCYHLT3:1:2146:50248:15716'

        
         

    def test_DAMID(self):

        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        index_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/'))
        barcode_parser = BarcodeParser(barcode_folder, lazyLoad='*',)
        index_parser = BarcodeParser(index_folder, lazyLoad='*')

        # First single cell format without overhang:
        # EG: DamID2_BC_001   3-TGCA-3-TATG
        
        first_barcode, second_barcode = 'TGCA', 'TATG'
        first_umi = 'ACT'
        second_umi='CTC'
        read_contents = 'rrrrrr'
        r1,r2 = self.construct_read_pair(f'{first_umi}{first_barcode}{second_umi}{second_barcode}',read_contents)
            
            
        damid_demux = DamID2_SCA(barcodeFileParser=barcode_parser,
                                    second_barcode_len=4,
                                    indexFileParser=index_parser,
                                    barcode_alias='DamID2_scattered_8bp',
                                        )

        demultiplexed_record = damid_demux.demultiplex([r1,r2])
        self.assertEqual( demultiplexed_record[0].tags['BC'], first_barcode+second_barcode)
        self.assertEqual( demultiplexed_record[0].tags['bi'], 1)
        self.assertEqual( demultiplexed_record[0].tags['RX'], first_umi+second_umi)
        self.assertEqual( demultiplexed_record[0].sequence, read_contents)
        
        combined_demux = DamID2andT_SCA(
            barcodeFileParser=barcode_parser,
            indexFileParser=index_parser)

        # The internal DamID demux of the combined protocol should return the same result
        demultiplexed_record = combined_demux.damid_demux.demultiplex([r1,r2])
        self.assertEqual( demultiplexed_record[0].tags['BC'], first_barcode+second_barcode)
        self.assertEqual( demultiplexed_record[0].tags['bi'], 1)
        self.assertEqual( demultiplexed_record[0].tags['RX'], first_umi+second_umi)
        self.assertEqual( demultiplexed_record[0].sequence, read_contents)
        
        
        demultiplexed_record = combined_demux.demultiplex([r1,r2])
        self.assertEqual( demultiplexed_record[0].tags['BC'], first_barcode+second_barcode)
        self.assertEqual( demultiplexed_record[0].tags['bi'], 1)
        self.assertEqual( demultiplexed_record[0].tags['dt'], 'DamID')
        self.assertEqual( demultiplexed_record[0].tags['RX'], first_umi+second_umi)
        self.assertEqual( demultiplexed_record[0].sequence, read_contents)


    def test_3DEC_UmiBarcodeDemuxMethod_matching_barcode(self):

        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        
        barcode_parser = BarcodeParser(barcode_folder,lazyLoad='*')


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

    def test_sra_header(self):

        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        barcode_parser = BarcodeParser(barcode_folder,lazyLoad='*')

        r1 = FastqRecord(
          '@SRR21016692.1 1/1',
          'ATCACACACTATAGTCATTCAGGAGCAGGTTCTTCAGGTTCCCTGTAGTTGTGTGGTTTTGAGTGAGTTTTTTAAT',
          '+',
          'AAAAA#EEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEE/EEEEEEEEEEEEEEEEEE'
        )
        r2 = FastqRecord(
          '@SRR21016692.1 1/2',
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



    def test_custom_header(self):

        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        barcode_parser = BarcodeParser(barcode_folder,lazyLoad='*')

        r1 = FastqRecord(
          '@CUSTOM_HEADER-WHICH_is_NOT-in_A_specific:format 1/1',
          'ATCACACACTATAGTCATTCAGGAGCAGGTTCTTCAGGTTCCCTGTAGTTGTGTGGTTTTGAGTGAGTTTTTTAAT',
          '+',
          'AAAAA#EEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEE/EEEEEEEEEEEEEEEEEE'
        )
        r2 = FastqRecord(
          '@CUSTOM_HEADER-WHICH_is_NOT-in_A_specific:format 1/2',
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
        self.assertEqual( demultiplexed_record[0].tags['oh'], 'CUSTOM_HEADER-WHICH_is_NOT-in_A_specific:format')

    def test_issue_290_header(self):

      barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
      barcode_parser = BarcodeParser(barcode_folder,lazyLoad='*')

      r1 = FastqRecord(
        '@LH00371:377:23JCYHLT3:1:1101:36621:1048 1:N:0:GACGAC',
        'CNGAAGGCTACTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTGGTGGGGGGTGTTTTTGTTTGAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
        '+',
        'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIII9III9I9I99II**III9I9I9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
      )
      r2 = FastqRecord(
        '@LH00371:377:23JCYHLT3:1:1101:36621:1048 1:N:0:GACGAC',
        'ANACGGCTAGGCCCTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTAGGGGGGGGGGTTGTGGTTTGAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
        '+',
        'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9I99I9999III9**9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
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
          indexFileAlias=None,
          random_primer_read=None,
          random_primer_length=6)

      demultiplexed_record = demux.demultiplex([r1,r2])
      assert 'oh' not in demultiplexed_record[0].tags

    
    def test_issue_290_header_provided_index(self):

      barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
      barcode_parser = BarcodeParser(barcode_folder,lazyLoad='*')

      index_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/'))
      index_parser = BarcodeParser(index_folder, lazyLoad='*', hammingDistanceExpansion=2)
      r1 = FastqRecord(
        '@LH00371:377:23JCYHLT3:1:1101:36621:1048 1:N:0:GACGAC',
        'CNGAAGGCTACTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTGGTGGGGGGTGTTTTTGTTTGAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
        '+',
        'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIII9III9I9I99II**III9I9I9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
      )
      r2 = FastqRecord(
        '@LH00371:377:23JCYHLT3:1:1101:36621:1048 1:N:0:GACGAC',
        'ANACGGCTAGGCCCTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTAGGGGGGGGGGTTGTGGTTTGAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
        '+',
        'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9I99I9999III9**9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
      )
      # Note that index GACGAC matches index number 41 in the index file illumina_RP_indices.bc
      demux = UmiBarcodeDemuxMethod(umiRead=0,
          umiStart=0,
          umiLength=3,
          barcodeRead=0,
          barcodeStart=3,
          barcodeLength=8,
          barcodeFileParser=barcode_parser,
          barcodeFileAlias='maya_384NLA',
          indexFileParser=index_parser,
          indexFileAlias='illumina_RP_indices',
          random_primer_read=None,
          random_primer_length=6)

      demultiplexed_record = demux.demultiplex([r1,r2])
      assert 'oh' not in demultiplexed_record[0].tags
      assert demultiplexed_record[0].tags['aa'] == 'GACGAC'
      assert demultiplexed_record[0].tags['aI'] == 41

        
    def test_issue_290_thruplex_header(self):
      # The same code as test_issue_290_header_provided_index, but then with a "normal" illumina_merged_ThruPlex48S_RP index
      barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
      barcode_parser = BarcodeParser(barcode_folder,lazyLoad='*')

      index_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/'))
      index_parser = BarcodeParser(index_folder, lazyLoad='*', hammingDistanceExpansion=2)
      r1 = FastqRecord(
        '@NB500901:415:HNVGWBGXK:4:11401:9263:1019 1:N:0:AGTTCC',
        'CNGAAGGCTACTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTGGTGGGGGGTGTTTTTGTTTGAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
        '+',
        'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIII9III9I9I99II**III9I9I9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
      )
      r2 = FastqRecord(
        '@NB500901:415:HNVGWBGXK:4:11401:9263:1019 1:N:0:AGTTCC',
        'ANACGGCTAGGCCCTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTAGGGGGGGGGGTTGTGGTTTGAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
        '+',
        'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9I99I9999III9**9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
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
      assert 'oh' not in demultiplexed_record[0].tags
      assert demultiplexed_record[0].tags['aa'] == 'AGTTCC'
      assert demultiplexed_record[0].tags['aI'] == 14



    

    def test_passing_indexFileAlias_to_baseDemux(self):
        """Test that indexFileAlias is correctly passed to baseDemux when handling rejected reads.
        
        This reproduces a bug where the DemultiplexingStrategyLoader creates baseDemux
        without passing the indexFileAlias, causing it to use the default instead of
        the user-specified one. When a read fails demultiplexing and gets written to
        the reject handle, index lookup fails if the index exists only in the specified
        file but not in the default.
        """
        
        barcode_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/'))
        barcode_parser = BarcodeParser(barcode_folder, lazyLoad='*')
        
        index_folder = str(importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/'))
        index_parser = BarcodeParser(index_folder, lazyLoad='*', hammingDistanceExpansion=2)
        
        # Create a read with GACGAC index (exists in illumina_RP_indices but NOT in default illumina_merged_ThruPlex48S_RP)
        # and a barcode that won't match (to trigger NonMultiplexable)
        r1 = FastqRecord(
            '@LH00371:377:23JCYHLT3:1:1101:36621:1048 1:N:0:GACGAC',
            'XXXNOTAREALBARCODEAAGCTACTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTGGTGGGGGGTGTTTTTGTTTGAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
            '+',
            'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIII9III9I9*I*99I*I**II*I*9I9I9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
        )
        r2 = FastqRecord(
            '@LH00371:377:23JCYHLT3:1:1101:36621:1048 2:N:0:GACGAC',
            'XXXNOTAREALBARCODEAAGCTAGGCCCTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTAGGGGGGGGGGTTGTGGTTTGAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
            '+',
            'I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII*I**III9I*99I*99**9*9III9**9*IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
        )
        
        # Create a temporary FASTQ file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(f"{r1.header}\n{r1.sequence}\n{r1.plus}\n{r1.qual}\n")
            temp_r1 = f.name
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(f"{r2.header}\n{r2.sequence}\n{r2.plus}\n{r2.qual}\n")
            temp_r2 = f.name
        
        try:
            # Create loader with illumina_RP_indices (not the default)
            loader = DemultiplexingStrategyLoader(
                barcodeParser=barcode_parser,
                indexParser=index_parser,
                indexFileAlias='illumina_RP_indices'  # NOT the default (illumina_merged_ThruPlex48S_RP)
            )
            
            # Create reject handle that can accept list writes (like real file handle)
            class ListRejectHandle:
                def __init__(self):
                    self.content = []
                def write(self, data):
                    if isinstance(data, list):
                        self.content.extend(data)
                    else:
                        self.content.append(data)
                def getvalue(self):
                    return '\n'.join(self.content)
            
            reject_output = ListRejectHandle()
            
            # This should not raise ValueError: Could not obtain index for GACGAC
            # The bug causes baseDemux to use default index instead of illumina_RP_indices
            processed, yields = loader.demultiplex(
                [temp_r1, temp_r2],
                rejectHandle=reject_output
            )
            
            # Verify that reads were rejected (because barcode doesn't match)
            self.assertEqual(processed, 1)  # One read pair processed
            
            # Verify the reject output was written without error
            # The header is now in tagged format with index info correctly populated
            reject_content = reject_output.getvalue()
            # Check that the index GACGAC was correctly resolved (aA and aI tags present)
            self.assertIn('aa:GACGAC', reject_content)
            self.assertIn('aA:GACGAC', reject_content)  # Corrected index
            self.assertIn('aI:41', reject_content)  # Index ID from illumina_RP_indices
            
        finally:
            os.unlink(temp_r1)
            os.unlink(temp_r2)


class TestTaggedFastqRoundTrip(unittest.TestCase):
    """Check if encode-decode round trip does not lose data"""

    def _parse(self, header):
        tr = TaggedRecord(TagDefinitions)
        tr.fromTaggedFastq(FastqRecord(header, 'ACGT', '+', 'IIII'))
        return tr.tags

    def test_roundtrip_tag_value_with_colons(self):
        # The fallback codec stores the original (non-illumina) header, colons included:
        oh = 'instrument:run:flowcell:lane:tile:x:y'
        tr = TaggedRecord(TagDefinitions)
        tr.tags.update({'Is': 'UNK', 'RN': '1', 'SM': 'lib_1', 'BC': 'ACACACTA',
                        'RX': 'ACGT', 'oh': oh})
        fastq = tr.asFastq(sequence='ACGT', dirAtt='+', baseQualities='IIII')
        parsed = self._parse(fastq.rstrip('\n').split('\n')[0])
        self.assertEqual(parsed['oh'], oh)
        # The other tags must not be affected:
        for tag in ('Is', 'RN', 'SM', 'BC', 'RX'):
            self.assertEqual(parsed[tag], tr.tags[tag])


    def test_roundtrip_colon_in_rejection_reason(self):
        reason = 'NoBarcodeMatch:ACGT'
        tr = TaggedRecord(TagDefinitions, reason=reason)
        fastq = tr.asFastq(sequence='ACGT', dirAtt='+', baseQualities='IIII')
        parsed = self._parse(fastq.rstrip('\n').split('\n')[0])
        self.assertEqual(parsed['RR'], reason)

    def test_segment_without_colon_is_skipped(self):
        # Legacy / malformed tail segments must not crash the parser:
        parsed = self._parse('@BC:ACACACTA;trailingjunk')
        self.assertEqual(parsed['BC'], 'ACACACTA')
        self.assertEqual(len(parsed), 1)

    def test_unregistered_key_does_not_crash(self):
        # Reject files written by the strategy loader carry the unregistered 'Rr' key:
        parsed = self._parse('@RR:not_a_reason;Rr:some error')
        self.assertEqual(parsed['RR'], 'not_a_reason')
        self.assertIn('Rr', parsed)

    def test_phred_tag_still_decoded(self):
        # QT is a registered phred tag: stored header-safe encoded by addTagByTag,
        # must be decoded back to ascii-phred by fromTaggedFastq:
        ascii_phred = 'IIIII'
        tr = TaggedRecord(TagDefinitions)
        tr.addTagByTag('QT', ascii_phred, isPhred=True)
        self.assertNotEqual(tr.tags['QT'], ascii_phred)  # encoded on write
        fastq = tr.asFastq(sequence='ACGT', dirAtt='+', baseQualities='IIII')
        parsed = self._parse(fastq.rstrip('\n').split('\n')[0])
        self.assertEqual(parsed['QT'], ascii_phred)

    def test_demultiplexed_fastq_conversion_tool(self):
        # End-to-end: demultiplexedFastqConversion.py must handle a file whose
        # headers carry an oh tag with colons (used to crash with ValueError):
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        script = os.path.join(repo_root, 'singlecellmultiomics',
                              'modularDemultiplexer', 'demultiplexedFastqConversion.py')
        tr = TaggedRecord(TagDefinitions)
        tr.tags.update({'BC': 'ACACACTA', 'QT': 'IIIIIIII',
                        'RX': 'ACGT', 'RQ': 'IIII',
                        'oh': 'instrument:run:flowcell:lane:tile:x:y'})
        with tempfile.NamedTemporaryFile('w', suffix='.fastq', delete=False) as handle:
            handle.write(tr.asFastq(sequence='NNNN', dirAtt='+', baseQualities='IIII'))
            temp_path = handle.name
        env = dict(os.environ, PYTHONPATH=repo_root)
        try:
            run = subprocess.run(
                [sys.executable, script,
                 '-seqfeatures', 'BC,RX', '-qualityfeatures', 'QT,RQ',
                 '-format', 'fq', temp_path],
                capture_output=True, text=True, env=env)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn('ACACACTAACGT', run.stdout)

            # -format is required now: unset must fail instead of silently emitting nothing
            run = subprocess.run(
                [sys.executable, script,
                 '-seqfeatures', 'BC,RX', '-qualityfeatures', 'QT,RQ', temp_path],
                capture_output=True, text=True, env=env)
            self.assertNotEqual(run.returncode, 0)

            # A missing tag must raise an error naming the tag:
            run = subprocess.run(
                [sys.executable, script,
                 '-seqfeatures', 'BC,NO', '-qualityfeatures', 'QT,RQ',
                 '-format', 'fq', temp_path],
                capture_output=True, text=True, env=env)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn('NO', run.stderr)
        finally:
            os.unlink(temp_path)


if __name__ == '__main__':
    unittest.main()
