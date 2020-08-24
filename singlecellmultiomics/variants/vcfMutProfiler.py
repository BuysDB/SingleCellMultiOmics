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
from singlecellmultiomics.variants.substitutions import conversion_dict, substitution_plot
import singlecellmultiomics


conversions_single_nuc = ["CA", "CG", "CT", "TA", "TC", "TG"]

def _get_conversions_per_cell(detected_variants_path, reference, known, prefix=None, **pysam_kwargs ):
    sc_detected = defaultdict(conversion_dict) # Cell -> patterncounts

    with pysam.VariantFile(detected_variants_path, **pysam_kwargs) as detected_vcf:
        try:
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

                        if prefix is not None:
                            sc_detected[(prefix,sample)][context, allele] += 1
                        else:
                            sc_detected[sample][context, allele] += 1
        except Exception as e:
            if not 'unable to parse' in str(e) and not 'has no len()' in str(e):
                raise
            else:
                return sc_detected

    return sc_detected


def get_conversions_per_cell( detected_variants_path, known_variants_path, reference_path, **kwargs ):

    reference = pysam.FastaFile(reference_path)

    known = set()
    if known_variants_path is not None:
        with pysam.VariantFile(known_variants_path) as known_vcf:
            for record in known_vcf.fetch():
                known.add( (record.chrom, record.pos ) )

    try:
        sc_detected = _get_conversions_per_cell(detected_variants_path, reference, known, **kwargs )
    except Exception as e:
        if str(e)=='no BGZF EOF marker; file may be truncated':
            print('no BGZF EOF marker; file may be truncated')
            sc_detected = _get_conversions_per_cell(detected_variants_path, reference, known, ignore_truncation=True, **kwargs )
        else:
            raise

    return sc_detected



if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Export and plot copy number profiles
    """)
    argparser.add_argument('detected_vcf',type=str, nargs='+')
    argparser.add_argument('-known',type=str)
    argparser.add_argument('-ref', help='path to reference fasta', type=str, required=True)
    argparser.add_argument('-rawmat', type=str, help='Path to write conversions table to (.csv / .pickle.gz)')
    argparser.add_argument('-heatmap', type=str, help='Path to write heatmap to (.png/.svg)')

    argparser.add_argument('-mut_count_threshold', default=0, type=int)

    argparser.add_argument('-bars', type=str, help='Path to write bargraphs to (.pdf)')

    args = argparser.parse_args()

    assert  args.rawmat is not None or args.heatmap is not None or args.bars is not None

    # Prune cells with very low counts:
    convs=None

    for vcf in args.detected_vcf:

        conversions = get_conversions_per_cell(
                                detected_variants_path = vcf,
                                known_variants_path = args.known,
                                reference_path=args.ref, prefix=(None if len(args.detected_vcf)==1 else vcf.split('/')[-1].replace('vcf','').replace('.gz','')) )
        if convs is None:
            convs = conversions
        else:
            convs.update(conversions)


    df = pd.DataFrame(convs)
    df = df.loc[:,df.sum()>=args.mut_count_threshold]

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

    if args.bars is not None:

        if args.mut_count_threshold>0:
            print('mut_count_threshold is ignored for -bars')

        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(args.bars,
                      metadata={'Creator': f'SingleCellMultiOmics {singlecellmultiomics.__version__}', 'Author': 'SCMO',
                                'Title': 'Mutation profiles'}) as pdf:



            for sample, pattern_counts in convs.items():
                print(sample)
                fig, ax = plt.subplots(figsize=(10,4))

                substitution_plot(pattern_counts, fig=fig, ax=ax, ylabel='# mutations')

                plt.suptitle(f'{sample}')
                pdf.savefig(fig)
                plt.close()




    print('All done')
