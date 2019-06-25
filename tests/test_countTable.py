#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from types import SimpleNamespace


import singlecellmultiomics.bamProcessing.bamToCountTable


"""
These tests check if the tagger is working correctly
"""

class TestMoleculeIteration(unittest.TestCase):

    def test_total_read_counting(self):
        df = singlecellmultiomics.bamProcessing.bamToCountTable.create_count_table(
            SimpleNamespace(
                alignmentfiles=['./data/mini_nla_test.bam'],
                o=None,
                bin=None,
                sliding=None,
                bedfile=None,
                showtags=False,
                featureTags=None,
                joinedFeatureTags='reference_name',
                sampleTags='SM',
                minMQ=0,
                filterXA=False,
                dedup=False,
                divideMultimapping=False,
                doNotDivideFragments=True,
                 noNames=False) , return_df=True)
        # !samtools idxstats ./singlecellmultiomics/data/mini_nla_test.bam | head -n 1 | cut -f 3
        self.assertEqual(df.loc['chr1'].sum(),563)

    def test_total_molecule_counting(self):
        df = singlecellmultiomics.bamProcessing.bamToCountTable.create_count_table(
            SimpleNamespace(
                alignmentfiles=['./data/mini_nla_test.bam'],
                o=None,
                bin=None,
                sliding=None,
                bedfile=None,
                showtags=False,
                featureTags=None,
                joinedFeatureTags='reference_name',
                sampleTags='SM',
                minMQ=0,
                filterXA=False,
                dedup=True,
                divideMultimapping=False,
                doNotDivideFragments=True,
                 noNames=False) , return_df=True)
        # !samtools view ./singlecellmultiomics/data/mini_nla_test.bam | grep 'RC:i:1' | wc -l
        self.assertEqual(df.loc['chr1'].sum(),383)


if __name__ == '__main__':
    unittest.main()
