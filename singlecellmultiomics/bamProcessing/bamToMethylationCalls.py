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
from singlecellmultiomics.methylation import MethylationCountMatrix


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


def get_methylation_count_matrix(bam_path,
                                 bin_size: int,
                                 bp_per_job: int,
                                 min_samples: int = None,
                                 min_variance: int = None,
                                 skip_contigs: set = None,
                                 known_variants: str = None,
                                 maxtime: int = None,
                                 threads: int = None):

    commands = generate_commands(
        alignments_path=bam_path,
        bin_size=bin_size,
        key_tags=None,
        max_fragment_size=0,
        dedup=True,
        bins_per_job=int(bp_per_job / bin_size), min_mq=40,
        kwargs={'known': known_variants,
                'maxtime': maxtime,
                }, skip_contigs=skip_contigs
    )


    count_mat = MethylationCountMatrix()
    if threads==1:
        for command in commands:
            result = count_methylation_binned(command)
            result.prune(min_samples=min_samples, min_variance=min_variance)
            count_mat.update( result )
    else:
        with multiprocessing.Pool(threads) as workers:

            for result in workers.imap_unordered(count_methylation_binned, commands):
                result.prune(min_samples=min_samples, min_variance=min_variance)
                count_mat.update(result)

    return count_mat



if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Export and plot copy number profiles
    """)
    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('-bin_size', default=500, type=int)
    argparser.add_argument('-bp_per_job', default=500_000, type=int)
    argparser.add_argument('-known_variants', help='VCF file with known variants, will be not taken into account as methylated/unmethylated', type=str)
    argparser.add_argument('-threads', default=None, type=int, help='Amount of threads to use, None to use the amount of available threads')

    fi = argparser.add_argument_group("Filters")
    fi.add_argument('-min_variance', default=0.1, type=float)
    fi.add_argument('-min_min_samples', default=10, type=int)

    og = argparser.add_argument_group("Output")
    og.add_argument('-bed', type=str, help='Bed file to write methylation calls to')
    og.add_argument('-wig', type=str, help='WIG file to write mean methylation per bin to')
    og.add_argument('-counts', type=str, help='CSV or pickle file to write single cell methylation counts to')
    og.add_argument('-bismark_tabfile', type=str, help='Tabulated file to write to, contains: chr | start | end | unmethylated_counts | methylated_counts | beta_value')
    og.add_argument('-tabfile', type=str,
                    help='Tabulated file to write to, contains: chr | start | end | unmethylated_counts | methylated_counts | beta_value | variance | n_samples')
    og.add_argument('-ref', help='path to reference fasta', type=str, required=True)
    og.add_argument('-counts', type=str, help='CSV or pickle file to write single cell methylation counts to')
    og.add_argument('-distmat', type=str, help='CSV or pickle file to write single cell distance matrix to')
    args = argparser.parse_args()

    counts = get_methylation_count_matrix(bam_path = args.bamfile,
                                 bin_size = args.bin_size,
                                 bp_per_job = args.bp_per_job,
                                 min_samples = args.min_samples,
                                 min_variance = args.min_variance,
                                 known_variants = args.known_variants,
                                 threads=args.threads,
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

    if args.distmat is not None:


        get_dmat(matrix)

    if args.bismark_tabfile is not None:
        counts.get_bulk_frame().to_csv(args.bismark_tabfile, sep='\t')

    if args.tabfile is not None:
        if args.tabfile.endswith('.csv'):
            counts.get_bulk_frame().to_csv(args.tabfile, sep='\t')
        else:
            counts.get_bulk_frame().to_pickle(args.tabfile)
