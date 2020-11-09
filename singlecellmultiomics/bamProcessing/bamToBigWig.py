#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyBigWig
import pysam
from singlecellmultiomics.bamProcessing.bamBinCounts import get_binned_counts
import pandas as pd
import argparse


def bam_to_wig(bam_paths, write_path, bin_size, method='sum', verbose=False):
    if verbose:
        print('Counting ...')
    counts = get_binned_counts(bam_paths, bin_size=bin_size)

    if verbose:
        print('Writing ...')
    with pysam.AlignmentFile(bam_paths[0]) as alignments, pyBigWig.open(write_path,'w') as out:

        # Write header
        out.addHeader(list(zip(alignments.references, alignments.lengths)))
        values = counts.sum(1).sort_index()

        # Write values
        for contig in alignments.references:
            if contig not in values.index.get_level_values(0):
                continue
            v = values.loc[[contig]]
            out.addEntries(
                list(v.index.get_level_values(0).values), #Contig
                list(v.index.get_level_values(1).values), #Start
                ends= list((v.index.get_level_values(1)+bin_size) .values), #end
                values= list(v.values))


if __name__ == '__main__':



    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Bam file to bigwig')

    argparser.add_argument('alignmentfiles', type=str, nargs='+')
    argparser.add_argument('-o', type=str, required=True, help='Output path (.bw)')

    argparser.add_argument('-bin_size', type=int, required=True)

    args = argparser.parse_args()

    bam_to_wig(args.alignmentfiles, args.o, args.bin_size, verbose=True)
