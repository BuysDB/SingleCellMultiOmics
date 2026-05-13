#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from singlecellmultiomics.fastqProcessing.fastqIterator import FastqRecord
from singlecellmultiomics.barcodeFileParser.barcodeFileParser import BarcodeParser
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.CELSeq2 import CELSeq2_c8_u6_NH
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.scCHIC import SCCHIC_384w_c8_u3_cs2
from singlecellmultiomics.modularDemultiplexer.demultiplexModules.DamID import DamID2andT_SCA,DamID2_SCA
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import TaggedRecord, UmiBarcodeDemuxMethod, TagDefinitions
from singlecellmultiomics.utils import reverse_complement
import importlib.resources


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
          indexFileAlias='illumina_merged_ThruPlex48S_RP',
          random_primer_read=None,
          random_primer_length=6)

      demultiplexed_record = demux.demultiplex([r1,r2])
      assert 'oh' not in demultiplexed_record[0].tags


    

if __name__ == '__main__':
    unittest.main()
