#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import itertools
import pysam
import singlecellmultiomics.universalBamTagger.universalBamTagger as ut
import singlecellmultiomics.universalBamTagger.bamtagmultiome as tm

"""
These tests check if the tagger is working correctly
"""

class TestMultiomeTaggingNLA(unittest.TestCase):

    def test_write_to_read_grouped_sorted(self):
        write_path = './data/write_test_rg.bam'
        tm.run_multiome_tagging_cmd(f'./data/mini_nla_test.bam -method nla -o {write_path}'.split(' '))

        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'bamtagmultiome' in x.get('PN','')]) )


if __name__ == '__main__':
    unittest.main()
