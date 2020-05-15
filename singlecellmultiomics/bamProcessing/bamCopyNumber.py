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
from singlecellmultiomics.bamProcessing.bamBinCounts import count_fragments_binned, generate_commands, gc_correct_cn_frame, obtain_counts
import os
import argparse
from colorama import Fore, Style


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Export and plot copy number profiles
    """)
    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('-ref', help='path to reference fasta', type=str, required=True)
    argparser.add_argument('-bin_size', default=500_000, type=int)
    argparser.add_argument('-max_cp', default=5, type=int)
    argparser.add_argument('-threads', default=16, type=int)
    argparser.add_argument('-bins_per_job', default=5, type=int)
    argparser.add_argument('-pct_clip', default=99.999, type=float)
    argparser.add_argument('-min_mapping_qual', default=40, type=int)
    argparser.add_argument('-molecule_threshold', default=5_000, type=int)

    argparser.add_argument('-rawmatplot', type=str, help='Path to raw matrix, plot is not made when this path is not supplied ')
    argparser.add_argument('-gcmatplot', type=str, help='Path to gc corrected matrix, plot is not made when this path is not supplied ')
    argparser.add_argument('-histplot', type=str, help='Path to histogram ')

    argparser.add_argument('-rawmat', type=str)
    argparser.add_argument('-gcmat', type=str)

    argparser.add_argument('-norm_method', default='median', type=str)

    args = argparser.parse_args()

    alignments_path = args.bamfile
    bin_size = args.bin_size
    MAXCP = args.max_cp
    pct_clip = args.pct_clip
    bins_per_job = args.bins_per_job
    min_mapping_qual = args.min_mapping_qual
    threads = args.threads
    molecule_threshold = args.molecule_threshold

    histplot=args.histplot
    rawmatplot=args.rawmatplot
    gcmatplot=args.gcmatplot
    rawmat=args.rawmat
    gcmat=args.gcmat

    reference = pysam.FastaFile(args.ref)
    h=GenomicPlot(reference)
    contigs = GenomicPlot(reference).contigs

    print("Creating count matrix ... ", end="")
    commands = generate_commands(
                alignments_path=alignments_path,
                bin_size=bin_size,key_tags=None,
                bins_per_job=5,head=None,min_mq=min_mapping_qual)

    counts = obtain_counts(commands,
                            reference=reference,
                            threads=threads,
                            live_update=False,
                            show_n_cells=None,
                            update_interval=None )
    print(f"\rCreating count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

    if histplot is not None:
        print("Creating molecule histogram ... ",end="")
        df = pd.DataFrame(counts).T.fillna(0)
        fig, ax = plt.subplots()
        cell_sums = df.sum()
        cell_sums.name = 'Frequency'
        cell_sums.plot.hist(bins=50)
        ax.set_xlabel('# molecules')
        ax.set_xscale('log')
        ax.axvline(molecule_threshold, c='r', label='Threshold')
        plt.legend()
        plt.savefig(histplot)
        print(f"\rCreating molecule histogram [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")


    # Convert the count dictionary to a dataframe
    print("Filtering count matrix ... ", end="")
    df = pd.DataFrame(counts).T.fillna(0)
    # remove cells were the median is zero
    if args.norm_method=='median':
        try:
            shape_before_median_filter = df.shape
            df = df.T[df.median()>0].T
            shape_after_median_filter = df.shape
            print(shape_before_median_filter,shape_after_median_filter )
            # Remove rows with little counts
            df = df.T[df.sum()>molecule_threshold].T
            df = df / np.percentile(df,pct_clip,axis=0)
            df = np.clip(0,MAXCP,(df / df.median())*2)
            df = df.T
        except Exception as e:
            print(f"\rMedian normalisation [ {Fore.RED}FAIL{Style.RESET_ALL} ] ")
            args.norm_method = 'mean'

    if args.norm_method == 'mean':
        shape_before_median_filter = df.shape
        df = df.T[df.mean()>0].T
        shape_after_median_filter = df.shape
        # Remove rows with little counts
        df = df.T[df.sum()>molecule_threshold].T
        df = df / np.percentile(df,pct_clip,axis=0)
        df = np.clip(0,MAXCP,(df / df.mean())*2)
        df = df.T



    if df.shape[0]==0:
        print(f"\rRaw count matrix [ {Fore.RED}FAIL{Style.RESET_ALL} ] ")
        raise ValueError('Resulting count matrix is empty, review the filter settings')
    else:
        print(f"\rFiltering count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")
        print( f'{df.shape[0]} cells, and {df.shape[1]} bins remaining' )
    del counts

    if rawmat is not None:
        print("Exporting raw count matrix ... ", end="")
        if rawmat.endswith('.pickle.gz'):
            df.to_pickle(rawmat)
        else:
            df.to_csv(rawmat)

        print(f"\rExporting raw count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

    if rawmatplot is not None:
        print("Creating raw heatmap ...", end="")
        h.cn_heatmap(df, figsize=(15,15))
        plt.savefig(rawmatplot)
        print(f"\rCreating raw heatmap [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")
        plt.close('all')

    if gcmatplot is not None or gcmat is not None:
        print("Performing GC correction ...", end="")
        corrected_cells = gc_correct_cn_frame(df, reference, MAXCP, threads, norm_method=args.norm_method)
        print(f"\rPerforming GC correction [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

    if gcmatplot is not None:
        print("Creating heatmap ...", end="")
        h.cn_heatmap(corrected_cells,figsize=(15,15))
        plt.savefig(gcmatplot)
        plt.close('all')
        print(f"\rCreating heatmap [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

    if gcmat is not None:
        print("Exporting corrected count matrix ... ")
        if gcmat.endswith('.pickle.gz'):
            corrected_cells.to_pickle(gcmat)
        else:
            corrected_cells.to_csv(gcmat)
        print(f"\rExporting corrected count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")
