#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itertools
import pysam
import os
import singlecellmultiomics.universalBamTagger.universalBamTagger as ut
import singlecellmultiomics.universalBamTagger.bamtagmultiome as tm

def test_write_to_read_grouped_multi():
    write_path = './data/write_test_chic_rg.bam'
    tm.run_multiome_tagging_cmd(f'./data/chic_test_region.bam -method chic --multiprocess -contig 8 -o {write_path}'.split(' '))

    with pysam.AlignmentFile(write_path) as f:
        # Test program header:

        i = 0
        # Test if the file has reads.
        for read in f:
            if read.is_read1:
                i += 1
            # Test if the reads have read groups:

    os.remove(write_path)
    os.remove(write_path + '.bai')

test_write_to_read_grouped_multi()