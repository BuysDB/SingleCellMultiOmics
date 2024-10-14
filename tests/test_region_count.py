#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from singlecellmultiomics.bamProcessing.bamCountRegions import count_regions
import tempfile
import pysam
from collections import Counter,defaultdict
import pandas as pd
import os

class TestRegionCount(unittest.TestCase):


    def test_region_counter(self):
        regions = {
            'regA':('chr1', 100,200),
            'regB':('chr1', 400,500),
            'regC':('chr1', 450,500),
        }

        reads = {
            # read name, cell index, library, mate, contig, cut location, strand : should be assigned  to:
            ('readA' , 1, 'LIBA', 1, 'chr1', 100, '+'): 'regA',
            ('readB' , 1, 'LIBA', 1, 'chr1', 150, '+'): 'regA',
            ('readC' , 1, 'LIBA', 1, 'chr1', 199, '+'): 'regA',
            ('readD' , 1, 'LIBA', 1, 'chr1', 200, '+'): None,
            ('readE' , 1, 'LIBA', 1, 'chr1', 210, '+'): None,
            ('readE' , 1, 'LIBA', 1, 'chr1', 399, '+'): None,
            ('readF' , 1, 'LIBA', 1, 'chr1', 400, '+'): 'regB',
            ('readG' , 1, 'LIBA', 1, 'chr1', 450, '+'): ['regB','regC'],
            ('readH' , 1, 'LIBB', 1, 'chr1', 450, '+'): ['regB','regC']

        }

        # Create bed file
        temp_bed_file = tempfile.NamedTemporaryFile(mode='w')
        for reg, (contig,start,end) in regions.items():
            temp_bed_file.write(f'{contig}\t{start}\t{end}\t{reg}\n')
        temp_bed_file.flush()

        # Create BAM file
        refseq = 'TTAATCATGAAACCGTGGAGGCAAATCGGAGTGTAAGGCTTGACTTTTAAGGGTTGAGATTTTCGAGAGGGATTCCTACGTTGCGTAGGTTCATGGGGGG'*10

        temp_bam_path = './temp_counting.bam'
        with pysam.AlignmentFile(
            temp_bam_path,
            'wb',
            reference_names=['chr1'],
            reference_lengths=[500]) as bam:

            for (read_name, cell_index, library, mate, contig, cut_location, strand) in reads:
                rlen = 20
                read_A = pysam.AlignedSegment(bam.header)
                read_A.reference_name = contig
                read_A.reference_start = cut_location
                # Before last A is a bogus G>A conversion to test strandness:
                read_A.query_sequence = refseq[cut_location:cut_location+rlen]
                read_A.cigarstring = f'{len(read_A.query_sequence)}M'
                read_A.qual = 'A'*len(read_A.query_sequence)
                read_A.mapping_quality = 60
                read_A.query_name = read_name
                read_A.set_tag('SM', f'{library}_{cell_index}')
                read_A.set_tag('LY', library)
                read_A.set_tag('DS', cut_location)
                read_A.set_tag('bi', cell_index)
                read_A.is_read1 = (mate==1)
                read_A.is_read2 = (mate!=1)

                read_A.is_paired = False
                read_A.is_proper_pair = True

                bam.write(read_A)
        pysam.index(temp_bam_path)

        ## Call the function to test:
        df = count_regions(temp_bed_file.name, bam_paths=[temp_bam_path,])
        temp_bed_file.close()


        targets = defaultdict(Counter)
        for (read_name, cell_index, library, mate, contig, cut_location, strand), goals in reads.items():
            if goals is None:
                continue
            if type(goals)!=list:
                goals = [goals]
            for goal in goals:
                targets[f'{library}_{cell_index}'][goal] += 1
        dfob = pd.DataFrame(targets).fillna(0)

        for feature, row in dfob.iterrows():
            for sample,value in row.items():
                self.assertTrue( df.loc[sample][feature] == value )

        os.remove(temp_bam_path)
        os.remove(temp_bam_path+'.bai')

if __name__ == '__main__':
    unittest.main()
