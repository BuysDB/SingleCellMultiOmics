#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import itertools

from singlecellmultiomics.features import FeatureContainer

"""
These tests check if the feature container is working correctly
"""

class TestFeatureContainer(unittest.TestCase):


    def expect(self, result, desired, presenceTestOnly=False):
        if presenceTestOnly:
            self.assertTrue( any( r[2] in desired for r in result) )
        if desired is None:
            self.assertTrue(len(result)==0)
        elif type(desired) is list:
            self.assertTrue( len(result)==len(desired) and all( r[2] in desired for r in result))
        else:
            self.assertTrue( len(result)==1 and result[0][2]==desired )

    def test_add(self):
        f = FeatureContainer()
        f.addFeature('chrY', 1, 3, 'A','+','')
        f.addFeature('chrY', 5, 8, 'B','+','')
        f.sort()
        self.expect( f.findFeaturesAt('chrY',0,'+'), None)
        self.expect( f.findFeaturesAt('chrY',1,'+'), 'A')
        self.expect( f.findFeaturesAt('chrY',2,'+'), 'A')
        self.expect( f.findFeaturesAt('chrY',3,'+'), 'A')
        self.expect( f.findFeaturesAt('chrY',4,'+'), None)
        self.expect( f.findFeaturesAt('chrY',5,'+'), 'B')
        self.expect( f.findFeaturesAt('chrY',6,'+'), 'B')
        self.expect( f.findFeaturesAt('chrY',8,'+'), 'B')

    def test_nested(self):
        f = FeatureContainer()
        f.addFeature('chrX', 10, 1000, 'parentB','+','')
        f.addFeature('chrX', 500, 900, 'nestedB','+','')
        f.addFeature('chrX', 10000, 12000, 'C','+','')
        f.addFeature('chrX', 100000, 120000, 'D','+','')
        f.sort()
        self.expect( f.findFeaturesAt('chrX',9,'+'), None)
        self.expect( f.findFeaturesAt('chrX',12001,'+'), None)
        self.expect( f.findFeaturesAt('chrX',12000,'+'), 'C')
        self.expect( f.findFeaturesAt('chrX',120000,'+'), 'D')
        self.expect( f.findFeaturesAt('chrX',10,'+'), 'parentB')

    def test_len(self):
        f = FeatureContainer()
        f.addFeature('chr1',100, 200,  '1','+','A forward feature from 100 to 200 chr1')
        f.addFeature('chr1',110, 200, '2','-',  'A reverse feature from 110 to 200 chr1')
        f.addFeature('chr2',100, 200,  '3', '+', 'A forward feature from 100 to 200 chr2')
        f.addFeature('chr2',100, 110,  '4','+', 'A forward feature from 100 to 110 chr2')
        f.addFeature('chr2',100, 150,  '5', '-','A reverse feature from 100 to 150 chr2')
        self.assertEqual( len(f) , 5)

    def test_stranded(self):
        f = FeatureContainer()
        f.addFeature('chr1',100, 200,  '1','+','A forward feature from 100 to 200 chr1')
        f.addFeature('chr1',110, 200, '2','-',  'A reverse feature from 110 to 200 chr1')
        f.addFeature('chr2',100, 200,  '3', '+', 'A forward feature from 100 to 200 chr2')
        f.addFeature('chr2',100, 110,  '4','+', 'A forward feature from 100 to 110 chr2')
        f.addFeature('chr2',100, 150,  '5', '-','A reverse feature from 100 to 150 chr2')

        f.addFeature('chr3',100, 150,  '6', '-','feature 6')
        f.addFeature('chr3',200, 250,  '7', '-','feature 7')
        f.addFeature('chr3',200, 450,  '8', '-','feature 8')
        f.addFeature('chr3',10, 15,  '9', '-','feature 9')
        f.sort()

        #printFormatted("[BRIGHT]Test for reference presence:")
        result = f.getReferenceList()
        desired = ['chr1','chr2','chr3']
        #printFormatted("self.expecting %s, %s" % (desired, ("[BRIGHT][GREEN] SUCCES [RESET][DIM]%s\n" % result) if len(result)==len(desired) and all( r in desired for r in result) else '[RED]FAIL %s\n' % result ))

        #printFormatted("[BRIGHT]Test on leftmost start:")
        self.expect( f.findFeaturesAt('chr1',100,'+'), '1')

        #printFormatted("[BRIGHT]Test on rightmost end:")
        self.expect( f.findFeaturesAt('chr1',200,'+'), '1')

        #printFormatted("[BRIGHT]Test on random location within feature:")
        self.expect(  f.findFeaturesAt('chr1',120,'+'), '1')

        #printFormatted("[BRIGHT]Test on random location within feature:")
        self.expect( f.findFeaturesAt('chr2',120,'+'), '3')

        #printFormatted("[BRIGHT]Test on limit location of feature:")
        self.expect( f.findFeaturesAt('chr2',200,'+'), '3')

        #printFormatted("[BRIGHT]Test on non matching location (match available on other side, and one base left of coord):")
        self.expect( f.findFeaturesAt('chr2',151,'-'), None)


    def test_double_matching(self):
        f = FeatureContainer()
        f.addFeature('chr1',100, 200,  '1','+','A forward feature from 100 to 200 chr1')
        f.addFeature('chr1',110, 200, '2','-',  'A reverse feature from 110 to 200 chr1')
        f.addFeature('chr2',100, 200,  '3', '+', 'A forward feature from 100 to 200 chr2')
        f.addFeature('chr2',100, 110,  '4','+', 'A forward feature from 100 to 110 chr2')
        f.addFeature('chr2',100, 150,  '5', '-','A reverse feature from 100 to 150 chr2')

        f.addFeature('chr3',100, 150,  '6', '-','feature 6')
        f.addFeature('chr3',200, 250,  '7', '-','feature 7')
        f.addFeature('chr3',200, 450,  '8', '-','feature 8')
        f.addFeature('chr3',10, 15,  '9', '-','feature 9')
        f.sort()
        #printFormatted("[BRIGHT]Tests on double matching locations without strand spec")
        self.expect( f.findFeaturesAt('chr1',120), ['1','2'])
        self.expect( f.findFeaturesAt('chr2',105), ['3','4','5'])


    def test_ranges(self):
        f = FeatureContainer()
        f.addFeature('chr1',100, 200,  '1','+','A forward feature from 100 to 200 chr1')
        f.addFeature('chr1',110, 200, '2','-',  'A reverse feature from 110 to 200 chr1')
        f.addFeature('chr2',100, 200,  '3', '+', 'A forward feature from 100 to 200 chr2')
        f.addFeature('chr2',100, 110,  '4','+', 'A forward feature from 100 to 110 chr2')
        f.addFeature('chr2',100, 150,  '5', '-','A reverse feature from 100 to 150 chr2')

        f.addFeature('chr3',100, 150,  '6', '-','feature 6')
        f.addFeature('chr3',200, 250,  '7', '-','feature 7')
        f.addFeature('chr3',200, 450,  '8', '-','feature 8')
        f.addFeature('chr3',10, 15,  '9', '-','feature 9')
        f.sort()

        #printFormatted("[BRIGHT] ==== Range tests... ====" )
        #printFormatted("[BRIGHT]Test for matching start and end coordinates overlapping one feature")
        self.expect(  f.findFeaturesBetween('chr2',102, 200, '+' ), ['3','4'])

        #printFormatted("[BRIGHT]Test for matching all features on chromosome")
        self.expect(  f.findFeaturesBetween('chr2',0, 20000, None ), ['3','4','5'])

        #printFormatted("[BRIGHT]Test for matching all but one features on chromosome")
        self.expect(  f.findFeaturesBetween('chr3',151, 20000, None ), ['8','7'])

        #printFormatted("[BRIGHT]Test for matching all but one features on chromosome")
        self.expect(  f.findFeaturesBetween('chr3',0, 195, None ), ['6','9'])

        #printFormatted("[BRIGHT]Test for range finding non-existent feature")
        self.expect(  f.findFeaturesBetween('chr2', 2001, 2000, '+' ), None)


        #printFormatted("[BRIGHT]Test for finding non-existent feature LEFT NEXT to the point")
        self.expect(  f.findNearestLeftFeature('chr2', 50, '+' ), None)

        #printFormatted("[BRIGHT]Test for finding existent feature LEFT NEXT to the point")
        self.expect(  f.findNearestLeftFeature('chr2', 250, None ), '3')


        #printFormatted("[BRIGHT]Test for finding non-existent feature RIGHT NEXT to the point")
        self.expect(  f.findNearestRightFeature('chr2', 250, '+' ), None)

        #printFormatted("[BRIGHT]Test for finding existent feature RIGHT NEXT to the point")
        self.expect(  f.findNearestRightFeature('chr2', 0, '-' ), '5')


        #printFormatted("[BRIGHT]Test for finding closest feature")
        self.expect(  f.findNearestFeature('chr1', 0, None ), '1')


if __name__ == '__main__':
    unittest.main()
