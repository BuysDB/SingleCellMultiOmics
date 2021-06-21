#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import itertools
import pysam
import os
import singlecellmultiomics.universalBamTagger.universalBamTagger as ut
import singlecellmultiomics.universalBamTagger.bamtagmultiome as tm

"""
These tests check if the tagger is working correctly
"""


class TestMultiomeTaggingCHIC(unittest.TestCase):

    def test_write_to_read_grouped_sorted(self):
        write_path = './data/write_test_chic_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/chic_test_region.bam -method chic -o {write_path}'.split(' '))


        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            # Test if the file has reads.
            for read in f:
                if read.is_read1:
                    i+=1
                self.assertTrue(read.has_tag('RG'))
            self.assertEqual(i, 17)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')

    def test_write_to_read_grouped_multi(self):
        write_path = './data/write_test_chic_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/chic_test_region.bam -contig 8 -method chic --multiprocess -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            # Test program header:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            # Test if the file has reads.
            for read in f:
                if read.is_read1:
                    i+=1
                # Test if the reads have read groups:
                self.assertTrue(read.has_tag('RG'))
            self.assertEqual(i, 17)




        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')

class TestMultiomeTaggingNLA(unittest.TestCase):

    def test_write_to_read_grouped_sorted(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam --allow_cycle_shift -method nla -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            qc_failed_R1 = 0
            # Test if the file has reads.
            for read in f:
                if read.is_read1:
                    i+=1
                    if read.is_qcfail:
                        qc_failed_R1+=1
                self.assertTrue(read.has_tag('RG'))
            self.assertEqual(i, 293)
            self.assertEqual(qc_failed_R1, 10)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')

    def test_tag_no_cycle_shift(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam -method nla -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            qc_failed_R1 = 0
            # Test if the file has reads.
            for read in f:
                self.assertTrue(read.has_tag('RG'))
                if read.is_read1:
                    i+=1
                    if read.is_qcfail:
                        qc_failed_R1+=1
            self.assertEqual(i, 293)
            self.assertEqual(qc_failed_R1, 13)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')

    def test_write_to_read_grouped_sorted_no_rejects(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam --no_rejects --allow_cycle_shift -method nla -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            # Test if the file has reads.
            for read in f:
                self.assertTrue(read.has_tag('RG'))
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 283)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')

    def test_skip_contig(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam --no_rejects --allow_cycle_shift -method nla -skip_contig chr1,chrMT -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            # Test if the file has reads.
            for read in f:
                self.assertTrue(read.has_tag('RG'))
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 0)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')

    def test_skip_contig_invert(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam --no_rejects --allow_cycle_shift -method nla -skip_contig chr2,chr3 -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            # Test if the file has reads.
            for read in f:
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 283)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')


    def test_skip_contig_multi_process(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam --multiprocess --no_rejects --allow_cycle_shift -method nla -skip_contig chr1,chrMT -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            # Test if the file has reads.
            for read in f:
                self.assertTrue(read.has_tag('RG'))
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 0)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')

    def test_skip_contig_invert_multi_process(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam --multiprocess --no_rejects --allow_cycle_shift -method nla -skip_contig chr2,chr3 -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )

            i =0
            # Test if the file has reads.
            for read in f:
                self.assertTrue(read.has_tag('RG'))
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 283)

        self.assertTrue( os.path.exists(write_path) )
        os.remove(write_path)
        os.remove(write_path+'.bai')


if __name__ == '__main__':
    unittest.main()
