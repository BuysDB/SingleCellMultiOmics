#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam
import collections
import pandas as pd
import numpy as np
import os


def obtain_approximate_reference_cut_position(site, contig, alt_spans):
    #contig, cut_start, strand = molecule.get_cut_site()
    alt_contig, alt_start, alt_end = alt_spans[contig]
    return contig, site + alt_start


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Create methylation and copy number tables from bam file, ALT contig aware')

    argparser.add_argument('alignmentfile', type=str)
    argparser.add_argument('-head', type=int, default=None)
    argparser.add_argument('-bin_size', type=int, default=250_000)
    argparser.add_argument(
        '-min_mq',
        type=int,
        default=30,
        help='Minimum mapping quality')
    argparser.add_argument('--met', action='store_true')
    argparser.add_argument('--cnv', action='store_true')
    argparser.add_argument('-ref', type=str, required=True)

    args = argparser.parse_args()

    #chromosomes = [f'chr{x}' for x in range(1,23)] + ['chrX']

    output_alias = args.alignmentfile + f'.table.bs_{args.bin_size}'
    # Load alternative contigs if available:
    alt_path = f'{args.ref}.64.alt'
    alt_spans = None
    if os.path.exists(alt_path):
        print(f'Loading ALT data from {alt_path}')
        with pysam.AlignmentFile(alt_path) as alts:
            alt_spans = {}
            for alt in alts:
                alt_spans[alt.query_name] = (
                    alt.reference_name, alt.reference_start, alt.reference_end)

    # cell -> (allele, chromosome,   bin) -> umi_count
    fragment_abundance = collections.defaultdict(collections.Counter)
    # cell -> (context, chromosome, bin) -> methylation
    methylation_pos = collections.defaultdict(collections.Counter)
    # cell -> (context, chromosome, bin) -> methylation
    methylation_neg = collections.defaultdict(collections.Counter)
    # cell -> (allele, chromosome,   bin) -> raw_count
    fragment_read_abundance = collections.defaultdict(collections.Counter)

    with pysam.AlignmentFile(args.alignmentfile, threads=4) as alignments:

        min_mq = 30
        wrote = 0
        # for chromosome in chromosomes:
        #    for i,read in enumerate(alignments.fetch(chromosome)):
        for read in alignments:
            if read.is_duplicate:
                continue
            if read.has_tag('mp') and read.get_tag('mp') != 'unique':
                continue
            if read.mapping_quality < args.min_mq or not read.has_tag('DS'):
                continue
            if args.head is not None and wrote >= (args.head - 1):
                break

            if read.has_tag('RZ') and read.get_tag('RZ') != 'CATG':
                continue

            contig = read.reference_name
            site = int(read.get_tag('DS'))
            if alt_spans is not None and contig in alt_spans:
                contig, site = obtain_approximate_reference_cut_position(
                    site, contig, alt_spans)

            wrote += 1
            bin_i = int(site / args.bin_size)

            allele = 'unk'
            if read.has_tag('DA'):
                allele = read.get_tag('DA')

            if args.met:
                for context in 'xhz':
                    if read.has_tag(f's{context}'):
                        methylation_neg[read.get_tag('SM')][(
                            context, allele, contig, bin_i)] += int(read.get_tag(f's{context}'))
                    if read.has_tag(f's{context.upper()}'):
                        methylation_pos[read.get_tag('SM')][(
                            context, allele, contig, bin_i)] += int(read.get_tag(f's{context.upper()}'))

            fragment_abundance[read.get_tag(
                'SM')][(allele, contig, bin_i)] += 1
            fragment_read_abundance[read.get_tag('SM')][(
                allele, contig, bin_i)] += read.get_tag('af')

    if args.cnv:
        print('writing count table')
        pd.DataFrame(fragment_abundance).to_pickle(
            f'{output_alias}.CNV_umis.pickle.gz')

        pd.DataFrame(fragment_read_abundance).to_pickle(
            f'{output_alias}.CNV_reads.pickle.gz')

    if args.met:
        print('writing count table for methylation status')
        methylation_pos = pd.DataFrame(
            methylation_pos).T.fillna(0).sort_index()
        methylation_neg = pd.DataFrame(
            methylation_neg).T.fillna(0).sort_index()

        for context in 'xhz':
            methylation_pos[context].T.to_pickle(
                f'{output_alias}.{context.upper()}.pickle.gz')
            methylation_neg[context].T.to_pickle(
                f'{output_alias}.{context}.pickle.gz')

            (methylation_pos[context] / (methylation_neg[context] + methylation_pos[context])
             ).T.to_pickle(f'{output_alias}.{context}.beta.pickle.gz')
