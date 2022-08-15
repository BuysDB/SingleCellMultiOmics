#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import numpy as np
import pysam
import seaborn as sns
from multiprocessing import Pool
from glob import glob
from singlecellmultiomics.utils import pool_wrapper
from singlecellmultiomics.bamProcessing.bamFunctions import get_samples_from_bam
import gzip

def _count_regions_cell_index(bamfile, regions, max_idx=384):
    counts = np.zeros((max_idx,len(regions)))
    library = None
    with pysam.AlignmentFile(bamfile) as alns:
        for i,(contig,start,end,label) in enumerate(regions):
            for read in alns.fetch(contig,start,end):
                if not read.is_read1 or read.is_qcfail or read.is_duplicate:
                    continue
                ds = read.get_tag('DS')
                if not (ds>=start and ds<=end):
                    continue
                idx = int(read.get_tag('bi'))
                counts[idx-1, i] += 1
        if read is not None and read.has_tag('LY'):
            library = read.get_tag('LY')

    if library is None:
        library = bamfile.split('/')[-2]
    indices = [f'{library}_{i}' for i in range(1,max_idx+1)]
    return pd.DataFrame(counts,
                        index=indices,
                        columns=[ label for (contig,start,end,label) in regions ])


def _count_regions(bamfile, regions):
    sample_list = get_samples_from_bam(bamfile)
    smap = {s:i for i,s in enumerate(sample_list)}
    counts = np.zeros((len(smap),len(regions)))

    with pysam.AlignmentFile(bamfile) as alns:
        for i,(contig,start,end,label) in enumerate(regions):
            for read in alns.fetch(contig,start,end):
                if not read.is_read1:
                    continue
                ds = read.get_tag('DS')
                if not (ds>=start and ds<=end):
                    continue

                idx = smap[read.get_tag('SM')]
                counts[idx, i] += 1

    return pd.DataFrame(counts,
                        index=sample_list,
                        columns=[ label for (contig,start,end,label) in regions ])

def count_regions(bed_path: str, bam_paths: list,  n_threads=None) -> pd.DataFrame:

    cmds = []
    with (gzip.open(bed_path,'rt') if bed_path.endswith('gz') else open(bed_path)) as o:
        for line in o:
            contig,start,end,tssid = line.strip().split()[:4]
            start,end = int(start), int(end)
            cmds.append((contig,start,end,tssid))

    dfs = []
    with Pool(n_threads) as workers:
        for df in workers.imap( pool_wrapper, (
            (_count_regions,
                 { 'bamfile': lib,
                   'regions': cmds,
                  # 'max_idx':max_idx
                   } )
            for lib in bam_paths
        )):
            dfs.append(df)

    return pd.concat( dfs )


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Count regions, uses DS tag. Only R1 reads are counted, qc_fail and duplicate reads are skipped.')

    argparser.add_argument('alignmentfiles', type=str, nargs='+')
    argparser.add_argument('-regions', type=str, required=True, help='bed file with regions to count')
    argparser.add_argument('-o', type=str, required=True, help='Output path (.csv, .csv.gz, .pickle.gz)')
    argparser.add_argument('-s', type=str, required=True, help='Grouped sum output path (.csv, .csv.gz, .pickle.gz)')
    argparser.add_argument('-t', type=int,  help='Amount of threads. Uses all available if not supplied')
    #argparser.add_argument('-max_idx', type=int, default=384, help='Maximum barcode index (the amount of wells in the plate)')

    args = argparser.parse_args()
    df = count_regions(args.regions, args.alignmentfiles, n_threads=args.t)

    if args.o.endswith('pickle.gz') or args.o.endswith('.pickle'):
        df.to_pickle(args.o)
    else:
        df.to_csv(args.o)

    if args.s is not None:
        df.columns = pd.MultiIndex.from_tuples( [col.split('_',1) for col in df.columns] )

        if args.s.endswith('pickle.gz') or  args.s.endswith('.pickle'):
            df.groupby(level=0,axis=1).sum().to_pickle(args.s)
        else:
            df.groupby(level=0,axis=1).sum().to_csv(args.s)
