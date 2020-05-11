#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
import seaborn as sns
import pysam
import numpy as np
import multiprocessing
from datetime import datetime
from singlecellmultiomics.utils.plotting import GenomicPlot
from singlecellmultiomics.bamProcessing.bamBinCounts import count_fragments_binned, generate_commands, gc_correct_cn_frame, obtain_counts, count_methylation_binned
import os
import argparse
from singlecellmultiomics.bamProcessing import get_samples_from_bam
from colorama import Fore, Style
from singlecellmultiomics.utils import dataframe_to_wig


def prefilter(counts, cell_names, min_measurements, min_variance):
    if min_measurements>0:
        counts = counts.loc[:, (counts >= 0).sum() > min_measurements]
    if min_variance>0:
        return counts.loc[:, counts.var() >= min_variance].reindex(cell_names)
    else:
        return counts.reindex(cell_names)


def panda_and_prefilter(args):
    d, args = args # counts_dict ,(cell_names, min_measurements, min_variance)
    return prefilter(pd.DataFrame(d), *args)


def get_methylation_count_matrix(bam_path, bin_size, bp_per_job,  min_measurements, min_variance, skip_contigs=['MT'],
                                 known_variants: str=None, maxtime=None, threads=None, output_fmt=None ):

    commands = generate_commands(
        alignments_path=bam_path,
        bin_size=bin_size,
        key_tags=None,
        max_fragment_size=0,
        dedup=True,
        bins_per_job=int(bp_per_job / bin_size), min_mq=40,
        kwargs={'known': known_variants,
                'maxtime': maxtime,
                'output_fmt':output_fmt
                }, skip_contigs=skip_contigs
    )

    # We need the cell names a priori for efficient memory management
    cell_names = get_samples_from_bam(bam_path)

    chunks = []

    if threads==1:
        for command in commands:
            chunks.append( count_methylation_binned(command) )
    else:
        with multiprocessing.Pool(threads) as workers:

            for result in workers.imap_unordered(count_methylation_binned, commands):
                chunks.append(result)


    with multiprocessing.Pool(threads) as workers:
        filtered_chunks = list(workers.imap(panda_and_prefilter, ( (c,(cell_names, min_measurements, min_variance)) for c in chunks)))

    # Record columns:
    columns = []
    for chunk in filtered_chunks:
        columns += list(chunk.columns)

    # Merge using numpy, because pandas eats the memory like a panda
    counts_matrix = np.concatenate(filtered_chunks, axis=1)
    del filtered_chunks
    del chunks

    return pd.DataFrame(counts_matrix, columns=pd.MultiIndex.from_tuples(columns), index=cell_names)



if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Export and plot copy number profiles
    """)
    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('-bin_size', default=500, type=int)
    argparser.add_argument('-bp_per_job', default=500_000, type=int)
    argparser.add_argument('-known_variants', help='VCF file with known variants', type=str)
    argparser.add_argument('-threads', default=None, type=int, help='Amount of threads to use, None to use the amount of available threads')
    argparser.add_argument('-min_variance', default=0.1, type=float)
    argparser.add_argument('-min_measurements', default=10, type=int)

    argparser.add_argument('-bed', type=str, help='Bed file to write methylation calls to')
    argparser.add_argument('-wig', type=str, help='WIG file to write mean methylation per bin to')
    argparser.add_argument('-counts', type=str, help='CSV or pickle file to write single cell counts to')
    argparser.add_argument('-bismark_tabfile', type=str, help='Tabulated file to write to, contains: chr | start | end | unmethylated_counts | methylated_counts | beta_value')


    argparser.add_argument('-ref', help='path to reference fasta', type=str, required=True)
    args = argparser.parse_args()


    counts = get_methylation_count_matrix(bam_path = args.bamfile,
                                 bin_size = args.bin_size,
                                 bp_per_job = args.bp_per_job,
                                 min_measurements = args.min_measurements,
                                 min_variance = args.min_variance,
                                 known_variants = args.known_variants,
                                 threads=args.threads
                                          )

    if args.counts is not None:
        if args.counts.endswith('.csv'):
            counts.to_csv(args.counts)
        else:
            counts.to_pickle(args.counts)


    if args.wig is not None:
        # Calculate mean beta value and write to wig:
        dataframe_to_wig(counts.mean().to_frame(), args.wig, span=args.bin_size)

    if args.bed is not None:
        pass

    if args.bismark_tabfile is not None:
        pass
