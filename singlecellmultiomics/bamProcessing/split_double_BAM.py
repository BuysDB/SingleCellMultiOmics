#!/usr/bin/env python
import sys, argparse, datetime
import collections
import os
import singlecellmultiomics
import collections
import itertools
import numpy as np
import random
import pysam
import pysamiterators
import matplotlib.colors
from importlib import reload
import pandas as pd
from scipy.interpolate import interp1d
from more_itertools import chunked
from more_itertools import windowed
import numpy as np

from singlecellmultiomics.bamProcessing import sorted_bam_file
from singlecellmultiomics.bamProcessing.bamToCountTable import coordinate_to_bins


def main():
    parser = argparse.ArgumentParser(description='Dual signal unmixing, through a probablity matrix for each cell across bins probability a read is assigned to signal 1. Use this prob matrix with bam file to split a bam file into signal 1 and signal 2')
    parser.add_argument('-inbam', metavar='INFILE',
                        help='Input bam file')
    parser.add_argument('-inprobmat', metavar='INFILE',
                        help='Tab sep matrix file. Columns are cell names (first fcolumn is ""). Rows are genomic bins. Values are probability of reads in bin assigned to mark1.')
    parser.add_argument('-outdir', metavar='OUTDIR',
                        help='Output directory for bams. Full name to be specified in script')
    parser.add_argument('-mapq', metavar='INTEGER 0 to 60', default=0, type=int,
                        help='Minimum quality of read to be considered')
    parser.add_argument('-binsize', metavar='Genomic binsize', default=50000, type=int,
                        help='Binsize of genomic bins to consider (assumes row names are defined by nearest 50kb bins)')
    parser.add_argument('--interpolation', action='store_true',
                        help='Makes a linear interpolation of the bins in your probability matrix (no interpolation across chromosomes).')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default=None,
                        help='Write arguments to logfile')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c'))
        print('Command line inputs:')
        print(CMD_INPUTS)
        print('Argparse variables:')
        print(ARG_INPUTS)


    p = pd.read_csv(args.inprobmat, sep="\t", index_col=0)

    def parse_bin_name(binname):
        chrname, coords = binname.split(':')
        start, end = coords.split('-')
        return chrname, int(start), int(end)

    if not args.interpolation:
        prob = p

    if args.interpolation:

        def interpolate_prob_mat(p):

            new_rows = []
            for index, (binA_orign, binB_orign) in enumerate(windowed(p.index, 2)):

                binA = binA_orign #parse_bin_name(binA_orign)
                binB = binB_orign #parse_bin_name(binB_orign)

                if binA[0] != binB[0]:
                    continue

                if binA[2] > binB[1]:
                    raise ValueError('The input is not sorted')

                contig = binA[0]

                binSize = binA[2] - binA[1]

                new_rows.append(p.loc[binA_orign, :])

                start, end = binA[2], binB[1]

                for new_bin_start in range(binA[2], binB[1], binSize):

                    new_bin_end = new_bin_start + binSize
                    new_bin_centroid = new_bin_start + binSize*0.5

                    # for every cell do interpolation
                    dx = end-start
                    d = (new_bin_centroid-start)
                    dy = p.loc[binB_orign, :] - p.loc[binA_orign, :]

                    interpolated = (dy/dx)*d + p.loc[binA_orign, :]
                    interpolated.name = (contig, new_bin_start, new_bin_end)

                    new_rows.append(interpolated)

            prob = pd.DataFrame(new_rows)

            indexNames = [f'{chromosomes}:{starts}-{ends}' for chromosomes, starts, ends in prob.index]
            prob.index = indexNames

            return prob

        p.index = pd.MultiIndex.from_tuples([parse_bin_name(t) for t in p.index])
        p = p.sort_index(0)

        prob = interpolate_prob_mat(p)

        prob.to_csv(os.path.join(args.outdir, "probabilityMatrix_linearInterpolated.csv"), sep='\t')

    #==========End interpolation============================================

    prob.index = pd.MultiIndex.from_tuples([parse_bin_name(t.replace('chr', '')) for t in prob.index])
    prob.index.set_names(["chr", "start", "end"], inplace=True)

    bamFile = args.inbam
    wrote = 0

    infboth = os.path.join(args.outdir, "both.bam")
    infA = os.path.join(args.outdir, "splitted_A.bam")
    infB = os.path.join(args.outdir, "splitted_B.bam")

    with pysam.AlignmentFile(bamFile) as f:
        with sorted_bam_file(infboth, f) as both, sorted_bam_file(infA, origin_bam=f) as a, sorted_bam_file(infB, origin_bam=f) as b:
            for readId, (R1, R2) in enumerate(pysamiterators.MatePairIterator(f)):
                if R1.mapping_quality < args.mapq & R2.mapping_quality < args.mapq:
                    continue  # one of two reads should have sufficient MAPQ. Less stringent. Should be OK?

                if R1.is_duplicate:
                    continue

                bin_start, bin_end = coordinate_to_bins(R1.get_tag('DS'), args.binsize, args.binsize)[0]

	    # Obtain prob:
                bin_name = (R1.reference_name, bin_start, bin_end)
                if not bin_name in prob.index:
                    continue
                if R1.get_tag('SM') not in prob.columns:
                    continue
                p = prob.loc[bin_name, R1.get_tag('SM')]
                wrote += 1
                group = 'A' if np.random.random() <= p else 'B'
                R1.set_tag('Gr', group)
                R2.set_tag('Gr', group)
                if group == 'A':
                    a.write(R1)
                    a.write(R2)
                else:
                    b.write(R1)
                    b.write(R2)
                both.write(R1)
                both.write(R2)
    print("Number of reads written:" + str(wrote))


if __name__ == '__main__':
    main()
