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



if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Export and plot copy number profiles
    """)
    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('-bin_size', default=500, type=int)
    argparser.add_argument('-bp_per_job', default=500_000, type=int)
    argparser.add_argument('-known_variants', help='VCF file with known variants', type=str)
    argparser.add_argument('-threads', default=None, type=int)
    argparser.add_argument('-ref', help='path to reference fasta', type=str, required=True)
    args = argparser.parse_args()

    bin_size = args.bin_size
    bp_per_job = args.bp_per_job




    print("Creating count matrix ... ", end="")


    def get_methylation_count_matrix(bam_path, bin_size, bp_per_job,  min_measurements, min_variance, skip_contigs=['MT'] ):

        commands = generate_commands(
            alignments_path=bam_path,
            bin_size=bin_size,
            key_tags=None,
            max_fragment_size=0,
            dedup=True,
            bins_per_job=int(bp_per_job / bin_size), min_mq=40,
            kwargs={'known': args.known_variants,
                    'maxtime': 3
                    }, skip_contigs=skip_contigs
        )

        # We need the cell names,
        cell_names = get_samples_from_bam(bam_path)

        chunks = []
        with multiprocessing.Pool(args.threads) as workers:

            for i, result in enumerate(workers.imap_unordered(count_methylation_binned,
                                                              commands)):
                chunks.append(result)


        def prefilter(counts, min_measurements=15, min_variance=0):
            counts = counts.loc[:, (counts >= 0).sum() > min_measurements]
            return counts.loc[:, counts.var() >= min_variance].reindex(cell_names)


        def panda_and_prefilter(d):
            return prefilter(pd.DataFrame(d))


        with multiprocessing.Pool(args.threads) as workers:
            filtered_chunks = list(workers.imap(panda_and_prefilter, chunks))

        # Record columns:
        columns = []
        for chunk in filtered_chunks:
            columns += list(chunk.columns)

        # Merge using numpy, because pandas eats the memory like a panda
        counts_matrix = np.concatenate(filtered_chunks, axis=1)
        del filtered_chunks
        del chunks

        counts = pd.DataFrame(counts_matrix, columns=pd.MultiIndex.from_tuples(columns), index=cell_names)
        return counts