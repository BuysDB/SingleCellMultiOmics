#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyBigWig
import pysam
from singlecellmultiomics.bamProcessing.bamBinCounts import get_binned_counts
import pandas as pd
import argparse
from singlecellmultiomics.bamProcessing import get_contig_sizes
import numpy as np

def bam_to_wig(bam_paths, write_path, bin_size, method='sum', verbose=False, n_threads=None):
    if verbose:
        print('Counting ...')
    counts = get_binned_counts(bam_paths, bin_size=bin_size, n_threads=n_threads)

    if verbose:
        print('Writing ...')
    with pysam.AlignmentFile(bam_paths[0]) as alignments, pyBigWig.open(write_path,'w') as out:

        cs = get_contig_sizes(alignments)
        # Write header
        out.addHeader(list(zip(alignments.references, alignments.lengths)))
        values = counts.sum(1).sort_index()
        print(values)
        # Write values
        for contig in alignments.references:

            if contig not in values.index.get_level_values(0):
                continue
            print(f'Writing data for {contig}')
            v = values.loc[[contig]]
            out.addEntries(
                list(v.index.get_level_values(0).values), #Contig
                list(v.index.get_level_values(1).values), #Start

                ends= list( np.clip(  (v.index.get_level_values(1)+bin_size) .values, 0, cs[contig]-1) ) ,  #end
                values= np.array(v.values, dtype=np.float32))


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Bam file to bigwig. Counts based on DS tag by default, otherwise falls back to the start coordinate of R1. Does not qcfail flagged reads and duplicate flagged reads.')

    argparser.add_argument('alignmentfiles', type=str, nargs='+')
    argparser.add_argument('-o', type=str, required=True, help='Output path (.bw)')

    argparser.add_argument('-bin_size', type=int, required=True)

    args = argparser.parse_args()

    bam_to_wig(args.alignmentfiles, args.o, args.bin_size, verbose=True)
