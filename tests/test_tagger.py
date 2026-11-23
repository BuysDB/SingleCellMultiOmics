#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import pysam
import os
import singlecellmultiomics.universalBamTagger.bamtagmultiome as tm
from bam_test_helpers import write_misc_alt_bam, count_reads_by_name

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


class TestMiscAltContigs(unittest.TestCase):
    """ Tests for alt contig (MISC_ALT_CONTIGS_SCMO) mode, unmapped read handling,
    multiprocess contig scope and region argument validation.
    A test bam is created with one main chromosome, two alt contigs and one unmapped read.
    """

    @classmethod
    def setUpClass(cls):
        # One main chromosome, two alt contigs and one unmapped read, see bam_test_helpers:
        cls.bam_path = write_misc_alt_bam('./data/misc_alt_test.bam')

    @classmethod
    def tearDownClass(cls):
        for path in [cls.bam_path, f'{cls.bam_path}.bai']:
            if os.path.exists(path):
                os.remove(path)

    def run_tagger(self, extra):
        """ Run the tagger and return the amount of times each read name was written """
        write_path = './data/write_test_misc.bam'
        tm.run_multiome_tagging_cmd([self.bam_path, '-method', 'qflag'] + extra + ['-o', write_path])
        counts = count_reads_by_name(write_path, ['main1', 'main2', 'alt1', 'alt2', 'alt3', 'unmapped1'])
        for path in [write_path, f'{write_path}.bai', write_path.replace('.bam', '.status.txt')]:
            if os.path.exists(path):
                os.remove(path)
        return counts

    def test_misc_alt_contigs_written_once(self):
        # Every alt contig molecule should be written exactly once, and no main or unmapped reads:
        self.assertEqual(self.run_tagger(['-contig', 'MISC_ALT_CONTIGS_SCMO']),
                         {'main1': 0, 'main2': 0, 'alt1': 1, 'alt2': 1, 'alt3': 1, 'unmapped1': 0})

    def test_unmapped_contig_written_once(self):
        # -contig '*' selects the unmapped reads only:
        self.assertEqual(self.run_tagger(['-contig', '*']),
                         {'main1': 0, 'main2': 0, 'alt1': 0, 'alt2': 0, 'alt3': 0, 'unmapped1': 1})

    def test_specific_contig_no_unmapped_reads(self):
        # Unmapped reads belong to no contig and should not be emitted by contig specific runs:
        self.assertEqual(self.run_tagger(['-contig', 'chr1']),
                         {'main1': 1, 'main2': 1, 'alt1': 0, 'alt2': 0, 'alt3': 0, 'unmapped1': 0})

    def test_full_genome_all_reads_once(self):
        self.assertEqual(self.run_tagger([]),
                         {'main1': 1, 'main2': 1, 'alt1': 1, 'alt2': 1, 'alt3': 1, 'unmapped1': 1})

    def test_misc_alt_contigs_multiprocess(self):
        self.assertEqual(
            self.run_tagger(['-contig', 'MISC_ALT_CONTIGS_SCMO', '--multiprocess', '-tagthreads', '2']),
            {'main1': 0, 'main2': 0, 'alt1': 1, 'alt2': 1, 'alt3': 1, 'unmapped1': 0})

    def test_specific_contig_multiprocess(self):
        self.assertEqual(
            self.run_tagger(['-contig', 'chr1', '--multiprocess', '-tagthreads', '2']),
            {'main1': 1, 'main2': 1, 'alt1': 0, 'alt2': 0, 'alt3': 0, 'unmapped1': 0})

    def test_region_with_contig(self):
        # main1 is located at 1000 and main2 at 2000 on chr1:
        self.assertEqual(
            self.run_tagger(['-contig', 'chr1', '-region_start', '900', '-region_end', '1500']),
            {'main1': 1, 'main2': 0, 'alt1': 0, 'alt2': 0, 'alt3': 0, 'unmapped1': 0})

    def test_region_requires_contig(self):
        # A region requires a single contig, and both region bounds need to be supplied:
        for extra in [['-region_start', '900', '-region_end', '1500'],              # no contig
                      ['-contig', 'chr1', '-region_end', '1500'],                   # no start coordinate
                      ['-contig', 'MISC_ALT_CONTIGS_SCMO',                          # pooled contigs have no coordinate space
                       '-region_start', '900', '-region_end', '1500']]:
            write_path = './data/write_test_misc_region.bam'
            with self.assertRaises(ValueError):
                tm.run_multiome_tagging_cmd([self.bam_path, '-method', 'qflag'] + extra + ['-o', write_path])
            # The status file is written before the arguments are validated:
            status_path = write_path.replace('.bam', '.status.txt')
            if os.path.exists(status_path):
                os.remove(status_path)


if __name__ == '__main__':
    unittest.main()
