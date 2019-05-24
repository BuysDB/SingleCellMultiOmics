#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import itertools

import singlecellmultiomics.barcodeFileParser.barcodeFileParser as barcodeFileParser


"""
These tests check if the barcode file parser is working correctly
"""

class TestBarcodeParser(unittest.TestCase):

    def test_insert(self):

        b = barcodeFileParser.BarcodeParser()

        b.addBarcode(
            barcodeFileAlias = 'test',
            barcode='AAA',
            index='TEST1')


        index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance('AAA','test')
        self.assertEqual(index,'TEST1')
        self.assertEqual(barcode,'AAA')
        self.assertEqual(hd,0)

    def test_hamming_expansion(self):

        b = barcodeFileParser.BarcodeParser()
        b.addBarcode( barcodeFileAlias = 'test', barcode='AAA', index='TEST1')
        b.expand(1,'test')
        # Postive examples:
        for positive in ['AAT','CAA','AGA','TAA']:
            index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance(positive,'test')
            self.assertEqual(index,'TEST1')
            self.assertEqual(barcode,'AAA')
            self.assertEqual(hd,1)

        #Negative examples
        for negative in ['TAT','CCA','GGA','TCA']:
            index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance(negative,'test')
            self.assertNotEqual(barcode,'AAA')

    def test_hamming_expansion(self):

        b = barcodeFileParser.BarcodeParser()
        b.addBarcode(
            barcodeFileAlias = 'test',
            barcode='AAA',
            index='TEST1')
        b.addBarcode(
            barcodeFileAlias = 'test',
            barcode='AAT',
            index='TEST2')
        b.addBarcode(
            barcodeFileAlias = 'test',
            barcode='TTT',
            index='TEST3')

        b.expand(2,'test')
        index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance('CTC','test')
        self.assertEqual(index,'TEST3')
        self.assertEqual(barcode,'TTT')
        self.assertEqual(hd,2)


    def test_hamming_expansion_collisions(self):

        b = barcodeFileParser.BarcodeParser()
        b.addBarcode(
            barcodeFileAlias = 'test',
            barcode='AAA',
            index='TEST1')
        b.addBarcode(
            barcodeFileAlias = 'test',
            barcode='AAT',
            index='TEST2')
        b.addBarcode(
            barcodeFileAlias = 'test',
            barcode='TTT',
            index='TEST3')

        b.expand(2,'test')
        index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance('AAG','test')
        self.assertIsNone(index)
        index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance('AAG','test')
        self.assertIsNone(index)
        index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance('TCA','test')
        self.assertIsNone(index)


        index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance('AAA','test')
        self.assertEqual(index,'TEST1')
        self.assertEqual(barcode,'AAA')
        self.assertEqual(hd,0)

        index, barcode, hd = b.getIndexCorrectedBarcodeAndHammingDistance('TAA','test')
        self.assertEqual(index,'TEST1')
        self.assertEqual(barcode,'AAA')
        self.assertEqual(hd,1)


if __name__ == '__main__':
    unittest.main()
