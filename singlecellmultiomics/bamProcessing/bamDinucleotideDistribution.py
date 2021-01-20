#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pandas as pd
from collections import Counter
from singlecellmultiomics.utils.sequtils import reverse_complement
from singlecellmultiomics.bamProcessing.bamAnalyzeCutDistances import get_sc_cut_dictionary
from singlecellmultiomics.bamProcessing import get_reference_path_from_bam, get_contigs_with_reads
from singlecellmultiomics.utils import is_main_chromosome
import os
from singlecellmultiomics.features import FeatureContainer
from itertools import product
import numpy as np
import pysam
from multiprocessing import Pool

mpl.rcParams['figure.dpi'] = 300


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Investigate dinucleotide usage')
    argparser.add_argument(
        '-o',
        type=str,
        help="output folder",
        default='./dinuc_usage/')
    argparser.add_argument(
            '-cpgislands',
            type=str,
            help="path to cpg island tsv file",
            required=False)

    argparser.add_argument(
            '-offset',
            type=float,
            help="Offset margin for the beta values, increase to get a larger y-axis",
            default=0.05)


    argparser.add_argument('alignmentfiles', type=str, nargs='+')

    args = argparser.parse_args()
    bams = args.alignmentfiles


    reference_path = get_reference_path_from_bam(args.alignmentfiles[0])

    # We need to ignore cpg islands for these stats:
    if args.cpgislands is None:
        cpg_islands = None
    else:
        cpg_islands = FeatureContainer()
        for idx,row in pd.read_csv(args.cpgislands,delimiter='\t').iterrows():
            cpg_islands.addFeature(chromosome=row['chrom'].replace('chr',''), start=row.chromStart,end= row.chromEnd, name=row.name)

        cpg_islands.sort()

    sc_cut_dict_stranded = get_sc_cut_dictionary(bams,strand_specific=True,)

    extraction_radius = 1500
    nucs = 'ACGT'
    dinucs =list(product(nucs,nucs))
    dinuc_map = dict( (d,i) for i,d in enumerate(dinucs))

    #reference = CachedFasta(reference_handle)
    def go(contig):
        nucobs = np.zeros([extraction_radius*2 + 1,4])
        dinucobs = np.zeros([extraction_radius*2 + 1,len(dinucs)])

        with pysam.FastaFile(reference_path) as reference:
            cut_obs = Counter()
            for cell, cellobs in sc_cut_dict_stranded[contig].items():
                cut_obs += Counter( list(cellobs.keys()) )

            for k, total_seen  in cut_obs.most_common():
                if total_seen<3: # Only use more common cuts, these are more likely relaxed
                    continue
                # Extract the genomic region:
                strand = k[0]
                position = k[1]

                if cpg_islands is not None:
                    hits = cpg_islands.findFeaturesAt(contig, position)
                    if len(hits):
                        continue

                if position-extraction_radius < 0:
                    continue

                seq = reference.fetch(contig, position-extraction_radius, position+extraction_radius)
                if strand:
                    seq = reverse_complement(seq)

                if not 'TA' in seq[extraction_radius-1:extraction_radius+1]: #or not 'TA' in seq[extraction_radius+145:extraction_radius+148]:
                    continue


                prevbase = None
                for i,base in enumerate(seq):
                    try:
                        j = nucs.index(base)
                    except ValueError:
                        prevbase = None
                        continue
                    try:
                        nucobs[i,j] += 1

                        k = dinuc_map.get((prevbase,base))
                        if k is not None:
                            dinucobs[i,k] += 1

                        prevbase= base
                    except IndexError:
                        prevbase = None
                        continue


        return nucobs, dinucobs

    if not os.path.exists(args.o):
        os.makedirs(args.o)

    nucobs = np.zeros([extraction_radius*2 + 1,4])
    dinucobs = np.zeros([extraction_radius*2 + 1,len(dinucs)])

    contigs  = [contig for contig in
        get_contigs_with_reads(bam_path=bams[0])
        if contig != 'Y' and contig[0]!='J' and contig!='MT' and is_main_chromosome(contig)]

    with Pool(16) as workers:
        for rn, rd in workers.imap(go, contigs):
            nucobs += rn
            dinucobs  += rd

    nuc_obs_df =  pd.DataFrame( nucobs, columns= list(nucs ))

    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 300

    def create_dinuc_plot(s=-250, e=170):

        things_for_legend = []
        fig,gcax = plt.subplots(figsize=(6,4))
        atax = gcax.twinx()

        dinuc_df = pd.DataFrame(dinucobs, columns= list(dinucs), index=list(range(-extraction_radius, extraction_radius+1)) )
        dinuc_df = dinuc_df.rolling(window=3,center=True).mean()


        dinuc_df = (dinuc_df.T / dinuc_df.sum(1)).loc[:,s:e].T

        gc = (dinuc_df[('C','C')] + dinuc_df[('G','G')] + dinuc_df[('G','C')])
        gcax.xaxis.label.set_color('orange')
        gcax.tick_params(axis='y', colors='orange')

        things_for_legend += gcax.plot(gc.index, gc.values, label='CC+GG+GC',c='orange')

        at  = (dinuc_df[('A','A')] + dinuc_df[('T','T')] + dinuc_df[('T','A')])
        atax.xaxis.label.set_color('purple')
        atax.tick_params(axis='y', colors='purple')

        things_for_legend += atax.plot( at.index, at.values, label='AA+TT+TA',c='purple')

        extra = args.offset
        atax.set_ylim(at.mean()-extra, at.mean()+extra)
        gcax.set_ylim(gc.mean()-extra, gc.mean()+extra)

        for x in range(s,e,10):
            atax.axvline(x+5,c='grey',lw=0.5, alpha=0.5)

        gcax.set_ylabel('Fraction CC+GG+GC')
        atax.set_ylabel('Fraction AA+AT+TA')
        gcax.set_xlabel('< into uncovered, Distance from cut (bp), into molecule >',c='k')
        atax.legend(things_for_legend, [t.get_label() for t in things_for_legend])
        plt.title('dinucleotide abundance')
        plt.axvline(0,c='k')
        plt.axvline(147/2,c='k')
        plt.axvline(147,c='k')

    create_dinuc_plot()
    plt.savefig(args.o + 'dinuc_usage_narrow.png', bbox_inches='tight')

    create_dinuc_plot(-500,500)
    plt.savefig(args.o + 'dinuc_usage.png', bbox_inches='tight')

    create_dinuc_plot(-1500,1500)
    plt.savefig(args.o + 'dinuc_usage_wide.png', bbox_inches='tight')

    dinuc_df = pd.DataFrame(dinucobs, columns= list(dinucs), index=list(range(-extraction_radius, extraction_radius+1)) )
    dinuc_df = dinuc_df.rolling(window=3,center=True).mean()
    dinuc_df = (dinuc_df.T / dinuc_df.sum(1)).T
    dinuc_df.to_pickle(args.o + 'dinucs.pickle.gz' )
    dinuc_df.to_csv(args.o + 'dinucs.csv' )
