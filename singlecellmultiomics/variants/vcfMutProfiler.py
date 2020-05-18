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

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
import itertools
from collections import OrderedDict, defaultdict
import seaborn as sns
from singlecellmultiomics.utils.sequtils import reverse_complement


conversions_single_nuc = ["CA", "CG", "CT", "TA", "TC", "TG"]

def get_conversions_per_cell( detected_variants_path, known_variants_path, reference_path ):

    reference = pysam.FastaFile(reference_path)

    known = set()
    if known_variants_path is not None:
        with pysam.VariantFile(known_variants_path) as known_vcf:
            for record in known_vcf.fetch():
                known.add( (record.chrom, record.pos ) )


    sc_detected = defaultdict(conversion_dict) # Cell -> patterncounts

    with pysam.VariantFile(detected_variants_path) as detected_vcf:
        for record in detected_vcf:
            if len(record.ref)!=1:
                continue

            if (record.chrom, record.pos) in known:
                continue

            origin_context = reference.fetch(record.contig, record.pos-2,record.pos+1).upper()
            for sample in record.samples:
                for allele in record.samples[sample].alleles:
                    if allele is None:
                        continue
                    if len(allele)!=1:
                        continue
                    if allele==record.ref:
                        continue


                    if not (record.ref+allele  in conversions_single_nuc):
                        context = reverse_complement(origin_context)
                        allele = reverse_complement(allele)
                    else:
                        context = origin_context
                    #print(allele,record.ref, context)
                    sc_detected[sample][context, allele] += 1


    return sc_detected

def conversion_dict():
    pattern_counts = OrderedDict()
    for ref, to in conversions_single_nuc:
        for context in itertools.product('ACGT',repeat=2 ):
            pattern_counts[(f'{context[0]}{ref}{context[1]}', to)] = 0
    return pattern_counts

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Export and plot copy number profiles
    """)
    argparser.add_argument('detected_vcf',type=str)
    argparser.add_argument('-known',type=str)
    argparser.add_argument('-ref', help='path to reference fasta', type=str, required=True)
    argparser.add_argument('-rawmat', type=str, help='Path to write conversions table to (.csv / .pickle.gz)')
    argparser.add_argument('-heatmap', type=str, help='Path to write heatmap to (.png/.svg)')
    argparser.add_argument('-mut_count_threshold', default=80, type=int)
    args = argparser.parse_args()

    assert  args.rawmat is not None or args.heatmap is not None

    # Prune cells with very low counts:
    df = pd.DataFrame(get_conversions_per_cell(
                                detected_variants_path = args.detected_vcf,
                                known_variants_path = args.known,
                                reference_path=args.ref ))
    df = df.loc[:,df.sum()>args.mut_count_threshold]

    # Write counts to disk:
    if  args.rawmat is not None:
        if args.rawmat.endswith('.pickle.gz'):
            df.to_pickle(args.rawmat)
        else:
            df.to_csv(args.rawmat)

        # Create heatmap
    if args.heatmap is not None:
        g = sns.clustermap( df,z_score=1, vmax=4, cmap='viridis',  yticklabels=True,method='ward' )
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 4)
        plt.savefig(args.heatmap)
        print(f"\rCreating heatmap [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")
        plt.close('all')


    print('All done')
