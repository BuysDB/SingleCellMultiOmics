#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import itertools
import pysam
import os
from singlecellmultiomics.alleleTools import AlleleResolver
import pysam

"""
These tests check if the AlleleResolver is working correctly
"""

class TestAlleleResolver(unittest.TestCase):

    def test_vcf_reader(self):

        test_vcf_path = './data/origin.vcf'
        vcf_string = """##fileformat=VCFv4.0
##reference=example.fa
##contig=<ID=1,length=42>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B
1\t18\t.\tA\tT\t42\tPASS\tDP=4\tGT\t1/1\t1/1
1\t20\t.\tA\tT\t42\tPASS\tDP=4\tGT\t0/0\t1/1
1\t22\t.\tG\tA\t42\tPASS\tDP=4\tGT\t0/0\t1/1
1\t40\t.\tA\tC\t42\tPASS\tDP=4\tGT\t./.\t1/1
"""
        with open(test_vcf_path,'w') as f:
            f.write(vcf_string)

        ar = AlleleResolver(vcffile=test_vcf_path)

        # Uninformative site:
        self.assertIsNone( ar.getAllelesAt('1',17,'A') )

        # Here are no matching alleles:
        self.assertIsNone( ar.getAllelesAt('1',19,'C') )

        # Sample A matches
        self.assertEqual( ar.getAllelesAt('1',19,'A'), set(['SAMPLE_A']) )

        # Sample B matches
        self.assertEqual( ar.getAllelesAt('1',21,'A'), set(['SAMPLE_B'] ) )

        # monomorphic: Sample B matches, sample A does not have the site
        self.assertEqual( ar.getAllelesAt('1',39,'C'), set(['SAMPLE_B'] ) )


        try:
            os.remove(test_vcf_path)
        except Exception as e:
            raise


if __name__ == '__main__':
    unittest.main()
