#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import itertools
import gzip
import os
import tempfile

import pysam

from singlecellmultiomics.features import FeatureContainer
from singlecellmultiomics.features.exonGTFtoIntronGTF import exonGTF_to_intronGTF

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


    def make_read(self, reference_start, cigarstring):
        """Create an aligned segment with the given alignment."""
        header = pysam.AlignmentHeader.from_dict(
            {'HD': {'VN': '1.6'}, 'SQ': [{'LN': 100000, 'SN': 'chr1'}]})
        read = pysam.AlignedSegment(header)
        read.query_name = 'test_read'
        read.reference_id = 0
        read.reference_start = reference_start
        read.cigarstring = cigarstring
        query_length = sum(
            length for op, length in read.cigartuples if op in (0, 1))
        read.query_sequence = 'A' * query_length
        read.qual = 'I' * query_length
        read.mapping_quality = 60
        return read

    def test_findFeaturesAtPysamAlign_abutting_blocks(self):
        """A read block ending exactly at a feature's first base must NOT hit it.

        pysam blocks are half-open [start, end), the container stores closed
        intervals (loadGTF offset=-1), so querying with the raw block end
        counted a phantom base one past the alignment.
        """
        f = FeatureContainer()
        # closed interval covering bases 100..199, as loadGTF stores exon 100-200:
        f.addFeature('chr1', 100, 199, 'ABUT', '+', None)
        f.sort()

        # block (50, 100) covers bases 50..99 only -> zero overlap:
        self.assertEqual(f.findFeaturesAtPysamAlign(self.make_read(50, '50M')), set())
        # one base further, block (50, 101) covers base 100 -> genuine overlap:
        self.assertEqual(
            {h[2] for h in f.findFeaturesAtPysamAlign(self.make_read(50, '51M'))},
            {'ABUT'})
        # abutting on the right was already rejected correctly before:
        self.assertEqual(f.findFeaturesAtPysamAlign(self.make_read(200, '10M')), set())

    def test_findFeaturesAtPysamAlign_methods_agree(self):
        """Base-wise (method 0) and block-wise (method 1) queries must agree."""
        f = FeatureContainer()
        f.addFeature('chr1', 40, 60, 'INSIDE_BLOCK1', '+', None)         # inside first block
        f.addFeature('chr1', 89, 95, 'LAST_BASE_BLOCK1', '+', None)      # touches last base of block 1
        f.addFeature('chr1', 90, 109, 'IN_GAP_ABUTTING', '+', None)      # intron bases only
        f.addFeature('chr1', 170, 200, 'ABUT_AFTER_BLOCK2', '+', None)   # starts 1 past last base
        f.sort()

        # blocks (50, 90) and (110, 170): covers bases 50..89 and 110..169:
        read = self.make_read(50, '40M20N60M')
        method0 = f.findFeaturesAtPysamAlign(read, method=0)
        method1 = f.findFeaturesAtPysamAlign(read, method=1)
        self.assertEqual(method0, method1)
        self.assertEqual(
            {h[2] for h in method0}, {'INSIDE_BLOCK1', 'LAST_BASE_BLOCK1'})

    def test_loadBED_closed_intervals(self):
        """BED is half-open and must be stored as closed intervals (like GTF)."""
        with tempfile.TemporaryDirectory() as tmpDir:
            path = os.path.join(tmpDir, 'regions.bed')
            with open(path, 'w') as bed:
                # region [50, 100) covers bases 50..99:
                bed.write('chr1\t50\t100\tSIMPLE\n')
                # blocked line: blocks covering 1000..1009 and 1050..1069:
                bed.write('chr1\t1000\t1100\tBLOCKED\t0\t+\t1000\t1100\t0\t2\t10,20,\t0,50,\n')
                # second block has zero length and covers nothing:
                bed.write('chr1\t2000\t2100\tZEROLEN\t0\t+\t2000\t2100\t0\t2\t10,0,\t0,50,\n')
            f = FeatureContainer()
            f.loadBED(path)

            def names_at(coord):
                return {h[2] for h in f.findFeaturesAt('chr1', coord)}

            self.assertEqual(names_at(50), {'SIMPLE'})
            self.assertEqual(names_at(99), {'SIMPLE'})
            # base 100 is just past the end of [50, 100):
            self.assertEqual(names_at(100), set())

            self.assertEqual(names_at(1000), {'BLOCKED'})
            self.assertEqual(names_at(1009), {'BLOCKED'})
            self.assertEqual(names_at(1010), set())
            self.assertEqual(names_at(1050), {'BLOCKED'})
            self.assertEqual(names_at(1069), {'BLOCKED'})
            self.assertEqual(names_at(1070), set())

            # zero-length block: the first block still hits, the empty one does not:
            self.assertEqual(names_at(2005), {'ZEROLEN'})
            self.assertEqual(names_at(2050), set())


class TestExonGTFtoIntronGTF(unittest.TestCase):

    gtfLines = [
        'chr1\ttest\texon\t100\t199\t.\t+\t.\ttranscript_id "T1"; gene_id "G1";\n',
        'chr1\ttest\texon\t300\t399\t.\t+\t.\ttranscript_id "T1"; gene_id "G1";\n',
        'chr2\ttest\texon\t601\t699\t.\t+\t.\ttranscript_id "T2"; gene_id "G2";\n',
        'chr2\ttest\texon\t801\t899\t.\t+\t.\ttranscript_id "T2"; gene_id "G2";\n',
        'chr2\ttest\texon\t1001\t1099\t.\t+\t.\ttranscript_id "T2"; gene_id "G2";\n',
        'chr3\ttest\texon\t50\t99\t.\t-\t.\ttranscript_id "T3"; gene_id "G3";\n',
        'chr3\ttest\texon\t200\t299\t.\t-\t.\ttranscript_id "T3"; gene_id "G3";\n',
    ]

    expectedIntrons = {
        ('chr1', 200, 299, '+'),
        ('chr2', 700, 800, '+'),  
        ('chr2', 900, 1000, '+'),
        ('chr3', 100, 199, '-'),
    }

    def writeGtf(self, path, gz=False):
        with (gzip.open(path, 'wt') if gz else open(path, 'w')) as handle:
            handle.writelines(self.gtfLines)

    def getIntrons(self, path):
        introns = set()
        for line in exonGTF_to_intronGTF(path, 'transcript_id'):
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            self.assertEqual(fields[2], 'intron')
            introns.add((fields[0], int(fields[3]), int(fields[4]), fields[6]))
        return introns

    def test_introns_across_chromosome_changes(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            for gz in (False, True):
                path = os.path.join(tmpDir, f'exons.gtf{".gz" if gz else ""}')
                self.writeGtf(path, gz)
                self.assertEqual(
                    self.getIntrons(path),
                    self.expectedIntrons,
                    msg=f'chromosome={path} gz={gz}'
                )

    def test_intronic_meta_data_passthrough(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            path = os.path.join(tmpDir, 'exons.gtf')
            self.writeGtf(path)
            for line in exonGTF_to_intronGTF(path, 'transcript_id'):
                fields = line.rstrip('\n').split('\t')
                if (fields[0], int(fields[3]), int(fields[4])) == ('chr2', 700, 800):
                    self.assertIn('T2', fields[8])
                    break
            else:
                self.fail('intron chr2:700-800 was not emitted')


if __name__ == '__main__':
    unittest.main()
