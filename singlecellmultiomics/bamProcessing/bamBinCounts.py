#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pysam
import collections
import numpy as np
import os
import pickle
import gzip
import pandas as pd
import multiprocessing
from singlecellmultiomics.bamProcessing import get_contig_sizes, get_contig_size
from singlecellmultiomics.utils import is_main_chromosome
from statsmodels.nonparametric.smoothers_lowess import lowess
from datetime import datetime

def obtain_approximate_reference_cut_position(site, contig, alt_spans):
    #contig, cut_start, strand = molecule.get_cut_site()
    alt_contig, alt_start, alt_end = alt_spans[contig]
    return contig, site + alt_start


def obtain_counts(commands, reference, live_update=True, show_n_cells=4, update_interval = 3, threads=4):

    if live_update:
        from singlecellmultiomics.utils.plotting import GenomicPlot
        import matplotlib.pyplot as plt
        cell_plots=  {}
        for cell_index in range(show_n_cells):
            gplot = GenomicPlot(reference)
            fig = gplot.get_figure()
            fig.canvas.draw()
            cell_plots[cell_index] = {'plot':gplot, 'fig':fig}

        plt.pause(0.01)

    #import random
    #random.shuffle(commands)

    counts = {}


    prev=None

    top_cells = None

    start_time=datetime.now()

    update_method='partial_df'

    with multiprocessing.Pool(threads) as workers:

        for i,result in enumerate(workers.imap_unordered(count_fragments_binned,
                              commands )):
            counts.update(result)
            if live_update and update_method=='partial_df':
                if (datetime.now()-start_time).total_seconds()>2 and (prev is None or (datetime.now()-prev).total_seconds() >= update_interval):
                    if len(result)==0:
                        continue

                    df = pd.DataFrame(counts).T
                    if df.sum().sum()==0:
                        continue
                    prev = datetime.now()

                    if top_cells is None:
                        top_cells = df.sum().sort_values()[-show_n_cells:]

                    df = df[top_cells.index].fillna(0)
                    df = np.clip(0,2, df / np.percentile(df,99,axis=0))

                    for contig in [list(result.keys())[0][0]]:
                        x = np.array([ (stop+start)/2 for start,stop in df.loc[contig].index.values] )

                        for cell_index, (cell, row) in enumerate(df.loc[contig].T.iterrows()):
                            gplot = cell_plots[cell_index]['plot']
                            gplot.reset_axis(contig)
                            fig = cell_plots[cell_index]['fig']
                            fig.suptitle(cell)
                            gplot[contig].scatter(x, row.values,s=0.1,c='k')
                            fig.canvas.draw()
                    plt.pause(0.001)

    # Show final result
    if live_update:
        df = pd.DataFrame(counts).T
        df = df[top_cells.index].fillna(0)
        df = np.clip(0,2, df / np.percentile(df,99,axis=0))

        for contig in cell_plots[0]['plot'].contigs:
            x = np.array([ (stop+start)/2 for start,stop in df.loc[contig].index.values] )

            for cell_index, (cell, row) in enumerate(df.loc[contig].T.iterrows()):
                gplot = cell_plots[cell_index]['plot']
                gplot.reset_axis(contig)
                fig = cell_plots[cell_index]['fig']
                fig.suptitle(cell)
                gplot[contig].scatter(x, row.values,s=0.1,c='k')
                fig.canvas.draw()
        plt.pause(0.001)

    return counts

def read_counts(read,min_mq, dedup=True):
    if not read.is_read1 or read is None:
        return False
    if dedup and read.is_duplicate:
        return False
    if read.has_tag('mp') and read.get_tag('mp') != 'unique':
        return False
    if read.mapping_quality < min_mq or not read.has_tag('DS'):
        return False
    #if read.has_tag('RZ') and read.get_tag('RZ') != 'CATG':
    #    return False
    return True


def gc_correct(args):
    observations, gc_vector,MAXCP = args
    correction = lowess(observations, gc_vector)
    return np.clip(observations / np.interp( gc_vector, correction[:,0], correction[:,1] ) , 0,MAXCP)

def gc_correct_cn_frame(df, reference, MAXCP, threads):

    # Perform GC correction
    chrom_sizes= dict( zip(reference.references, reference.lengths))
    # Extract GC percentage from reference for the selected bin size:
    bins_to_gc = {}

    for contig,start,end in df.columns:
        k =  (contig,start,end)
        if not k in bins_to_gc:
            sequence = reference.fetch(contig, start,end ).upper()
            gc = sequence.count('G')+sequence.count('C')
            gcrat = (gc) / ((sequence.count('A')+sequence.count('T')+gc))
            bins_to_gc[ k ] = gcrat

    # Join the GC table with the count matrix
    gc_matched = df.T.join( pd.DataFrame({'gc':bins_to_gc}), how='left')['gc']

        # This performs GC correction for every cell using loess regression
    with multiprocessing.Pool(threads) as workers:
        keep_bins=df.columns
        gc_vector = gc_matched[keep_bins]

        corrected_cells = list( workers.imap(
            gc_correct, [(row,gc_vector.values,MAXCP) for cell,row in df.iterrows()] ))

    corrected_cells = pd.concat(corrected_cells,axis=1).T
    corrected_cells = ((corrected_cells.T/corrected_cells.median(1))*2).T

    return corrected_cells



def generate_jobs(alignments_path, bin_size = 1_000_000, bins_per_job = 10):
    for job_group in (((contig, start, start+bin_size*bins_per_job)
                 for start in range(0,length,bin_size*bins_per_job))
                for contig,length in
                get_contig_sizes(alignments_path).items()):
        yield from job_group

def generate_commands(alignments_path, bin_size = 1_000_000, bins_per_job = 10,alt_spans=None, min_mq=50,max_fragment_size=1000, head=None,key_tags=None,dedup=True):
    for i,(contig,start,end) in enumerate(generate_jobs(alignments_path=alignments_path,bin_size=bin_size,bins_per_job=bins_per_job)):
        yield (alignments_path, bin_size, max_fragment_size, \
                               contig, start, end, \
                               min_mq,alt_spans,key_tags,dedup)
        if head is not None and i>=(head-1):
            break

def count_fragments_binned(args):
    (alignments_path, bin_size, max_fragment_size, \
                           contig, start, end, \
                           min_mq,alt_spans, key_tags,dedup) = args


    counts = {} # Sample->(contig,bin_start,bin_end)->counts
    # Define which reads we want to count:


    p=0
    with pysam.AlignmentFile(alignments_path, threads=4) as alignments:
        # Obtain size of selected contig:
        contig_size = get_contig_size(alignments,contig)
        if contig_size is None:
            raise ValueError('Unknown contig')

        # Determine where we start looking for fragments:
        f_start = max(0, start-max_fragment_size)
        f_end = min( end+max_fragment_size, contig_size)

        for p,read in enumerate(alignments.fetch(contig=contig,start=f_start,
                                                 stop=f_end)):

            if not read_counts(read,min_mq=min_mq, dedup=dedup):
                continue

            # Extract the site
            site = int(read.get_tag('DS'))

            # Don't count sites outside the selected bounds
            if site<start or site>=end:
                continue

            sample = read.get_tag('SM')

            # Process alternative contig counts:
            if alt_spans is not None and contig in alt_spans:
                contig, site = obtain_approximate_reference_cut_position(
                    site, contig, alt_spans)

            # Obtain the bin index
            bin_i = int(site / bin_size)
            bin_start = bin_size*bin_i
            bin_end = min( bin_size*(bin_i+1), contig_size )

            # Add additional tag information: (For example the allele tag)
            if key_tags is not None:
                tag_values = [ (read.get_tag(tag) if read.has_tag(tag) else None) for tag in key_tags ]
                bin_id = (*tag_values, contig, bin_start, bin_end)
            else:
                bin_id = (contig, bin_start, bin_end)

            # Add a (single) count tot the dictionary:
            if not bin_id in counts:
                counts[bin_id] = {}

            if not sample in counts[bin_id]:
                counts[bin_id][sample] = 1
            else:
                counts[bin_id][sample] += 1

    return counts


def count_fragments_binned_wrap(args):

    (alignments_path, bin_size, max_fragment_size, \
                           contig, start, end, \
                           min_mq,alt_spans,dedup) = args

    tp = f'./TEMP_{contig}_{start}.pickle.gz'
    res = os.system(f'bamBinCounts.py {alignments_path} -o {tp} -start {start} -end {end} -contig {contig}')
    with gzip.open(tp,'rb') as pf:
        result = pickle.load(pf)
    os.remove(tp)
    return result

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      description='Extract counts from bam file')

    argparser.add_argument('alignmentfile', type=str)
    argparser.add_argument('-head', type=int, default=None)
    argparser.add_argument('-bin_size', type=int, default=250_000)
    argparser.add_argument('-t', type=int, default=16, help='Threads')
    argparser.add_argument('-j', type=int, default=5, help='Bins per worker')
    argparser.add_argument(
      '-min_mq',
      type=int,
      default=30,
      help='Minimum mapping quality')
    argparser.add_argument('-max_fragment_size', type=int, default=1000)
    argparser.add_argument('-o', type=str, required=True)
    args = argparser.parse_args()

    counts = {}

    with multiprocessing.Pool(args.t) as workers:

        job_gen = generate_commands(
                alignments_path=args.alignmentfile,
                bin_size=args.bin_size,
                bins_per_job=args.j,
            head=args.head)

        for job_total,_ in enumerate(job_gen):
            pass

        job_gen = generate_commands(
                alignments_path=args.alignmentfile,
                bin_size=args.bin_size,
                bins_per_job=args.j,
            head=args.head)

        for i,result in enumerate(workers.imap_unordered(count_fragments_binned,
                              job_gen )):
            counts.update(result)
            print(f'\r{i}/{job_total}',end='')

        if args.o.endswith('.pickle.gz'):
            pd.DataFrame(counts).to_pickle(args.o)
        elif args.o.endswith('.csv'):
            pd.DataFrame(counts).to_csv(args.o)
        else:
            with gzip.open(args.o,'wb') as out:
                pickle.dump(counts, out)
