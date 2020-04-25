#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import singlecellmultiomics
from singlecellmultiomics.molecule import MoleculeIterator, NlaIIIMolecule
from singlecellmultiomics.fragment import NlaIIIFragment
import pysam
import collections
import pysamiterators
import pandas as pd
import seaborn as sns
import numpy as np
import re
import sklearn.ensemble
import argparse


def bin_to_sort_value(chrom):
    chrom = chrom.replace('chr', '')
    try:
        if '_' in chrom:
            int_chrom = int(chrom.split('_')[0])
        elif chrom == 'X':
            int_chrom = 99
        else:
            int_chrom = int(chrom)
    except Exception as e:
        return ((999, chrom))
    return (int_chrom, chrom)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Correct count table for GC, site abundance and others""")
    argparser.add_argument('-reads', required=True)
    argparser.add_argument('-umis', required=True)
    argparser.add_argument('-no-change-regions', required=True)
    argparser.add_argument('-ref', required=True)
    argparser.add_argument('-corrected_umis_path', required=True)
    argparser.add_argument('-corrected_reads_path', required=True)
    argparser.add_argument('-bin_size', required=True, type=int)

    args = argparser.parse_args()
    bin_size = args.bin_size


    no_change_regions = args.no_change_regions.split(',')

    print("Reading data...")
    df_reads = pd.read_pickle(args.reads).sum(level=[1,2], axis=0).T
    df_umis = pd.read_pickle(args.umis).sum(level=[1,2], axis=0).T

    df_reads = df_reads.reindex(
        sorted(
            df_reads.columns,
            key=lambda a: (
                bin_to_sort_value(
                    a[0]),
                a[1])),
        axis=1)

    df_umis = df_umis.reindex(
        sorted(
            df_umis.columns,
            key=lambda a: (
                bin_to_sort_value(
                    a[0]),
                a[1])),
        axis=1)

    print(f"data size: {df_umis.shape}")
    reference = pysamiterators.CachedFasta(pysam.FastaFile(args.ref))

    # Calculate reference statistics
    print("Calculating reference statistics ... ")
    ref_stats = collections.defaultdict(dict)
    for chrom, bin_idx in df_reads:
        bin_seq = reference.fetch(
            chrom,
            bin_idx * bin_size,
            (1 + bin_idx) * bin_size).upper()
        base_obs = collections.Counter(bin_seq)

        ref_stats[(chrom, bin_idx)]['BASES'] = sum(base_obs.values())
        if ref_stats[(chrom, bin_idx)]['BASES'] > 0:
            ref_stats[(chrom, bin_idx)]['GC'] = (base_obs['G'] +
                                                 base_obs['C']) / ref_stats[(chrom, bin_idx)]['BASES']
        else:
            # add the mean here later
            ref_stats[(chrom, bin_idx)]['GC'] = np.nan
        ref_stats[(chrom, bin_idx)]['SITE_COUNT'] = bin_seq.count('CATG')

        frag_sizes = np.diff([m.start() for m in re.finditer('CATG', bin_seq)])
        ref_stats[(chrom, bin_idx)]['FRAG>40'] = np.sum(frag_sizes > 40)
        ref_stats[(chrom, bin_idx)]['FRAG>70'] = np.sum(frag_sizes > 70)
        ref_stats[(chrom, bin_idx)]['MEAN_FS'] = np.mean(frag_sizes)


    regressor = sklearn.ensemble.RandomForestRegressor(
        n_estimators=100, n_jobs=8)

    print("Training model ... ")
    rdf = pd.DataFrame(ref_stats)
    X = rdf[no_change_regions].T
    y = df_reads[no_change_regions].sum(0) * 2

    regressor.fit(X, y)
    # fill with other value @todo
    print(f"Performing correction on {len(X)} bins ")
    reduced = (df_umis / regressor.predict(rdf.T)).fillna(0)
    pd.DataFrame(reduced).to_pickle(args.corrected_umis_path)

    reduced = (df_umis / regressor.predict(rdf.T)).fillna(0)
    pd.DataFrame(reduced).to_pickle(args.corrected_reads_path)
