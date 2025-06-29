#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics
import importlib.resources
"""
These tests check if the Molecule module is working correctly
"""

class TestAuxData(unittest.TestCase):

    def test_barcode_files_presence(self):
        """Test if the barcode files can be accessed"""
        available_barcodes =  [p.name  for p in (importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/barcodes/')).iterdir() if p.is_file()]

        at_least_required = ['nla_bisulfite.bc',
         'lennart96NLA.bc',
         'illumina_RP_indices.bc',
         'celseq1.bc',
         'maya_mspj1.bc',
         'scartrace.bc',
         'maya_384NLA.bc',
         'scartraceBarcodes.bc',
         'lk_virus1.bc',
         'celseq2_noNla.bc',
         'celseq2.bc']

        for a in at_least_required:
            self.assertIn(a,available_barcodes)

    def test_index_files_presence(self):
        """Test if the index files can be accessed"""
        available_barcodes =  [p.name  for p in (importlib.resources.files('singlecellmultiomics').joinpath('modularDemultiplexer/indices/')).iterdir() if p.is_file()]

        at_least_required = ['illumina_TruSeq_indices.bc',
             'illumina_RP_indices.bc',
             'illumina_merged_iPCR_RP.bc',
             'illumina_ThruPlex48S_indices.bc',
             'illumina_merged_ThruPlex48S_RP.bc',
             'illumina_i7_indices.bc']

        for a in at_least_required:
            self.assertIn(a,available_barcodes)

if __name__ == '__main__':
    unittest.main()
