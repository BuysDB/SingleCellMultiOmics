#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam
import collections
import pandas as pd
import numpy as np

argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Create methylation and copy number tables from bam file')

argparser.add_argument('alignmentfile',  type=str)
argparser.add_argument('-head',  type=int, default=None)
argparser.add_argument('-bin_size',  type=int,default=250_000)
argparser.add_argument('-min_mq',  type=int,default=30, help='Minimum mapping quality')
argparser.add_argument('--met',  action='store_true')
argparser.add_argument('--cnv',  action='store_true')
args = argparser.parse_args()

chromosomes = [f'chr{x}' for x in range(1,23)] + ['chrX']

output_alias = args.alignmentfile + f'.table.bs_{args.bin_size}'

fragment_abundance = collections.defaultdict( collections.Counter) # cell -> (allele, chromosome,   bin) -> umi_count
methylation_pos = collections.defaultdict( collections.Counter) # cell -> (context, chromosome, bin) -> methylation
methylation_neg = collections.defaultdict( collections.Counter) # cell -> (context, chromosome, bin) -> methylation

with pysam.AlignmentFile(args.alignmentfile) as alignments:

    min_mq = 30
    wrote = 0
    for chromosome in chromosomes:
        for i,read in enumerate(alignments.fetch(chromosome)):
            if read.is_duplicate:
                continue
            if read.mapping_quality<args.min_mq or not read.has_tag('DS'):
                continue
            if args.head is not None and wrote>=(args.head-1):
                break

            wrote +=1
            bin_i =  int( read.get_tag('DS') / args.bin_size )

            allele='unk'
            if read.has_tag('DA'):
                allele = read .get_tag('DA')


            if args.met:
                for context in 'xhz':
                    if read.has_tag(f's{context}'):
                        methylation_neg[read.get_tag('SM')][(context, allele, read.reference_name, bin_i)]+=int(read.get_tag(f's{context}'))
                    if read.has_tag(f's{context.upper()}'):
                        methylation_pos[read.get_tag('SM')][(context, allele, read.reference_name, bin_i)]+=int(read.get_tag(f's{context.upper()}'))

            fragment_abundance[read.get_tag('SM')][(allele, read.reference_name, bin_i)]+=1


if args.cnv:
    print('writing count table')
    pd.DataFrame( fragment_abundance ).to_pickle(f'{output_alias}.CNV.pickle.gz')


if args.met:
    print('writing count table for methylation status')
    methylation_pos = pd.DataFrame( methylation_pos ).T.fillna(0).sort_index()
    methylation_neg = pd.DataFrame( methylation_neg ).T.fillna(0).sort_index()

    for context in 'xhz':
        methylation_pos[context].T.to_pickle(f'{output_alias}.{context.upper()}.pickle.gz')
        methylation_neg[context].T.to_pickle(f'{output_alias}.{context}.pickle.gz')

        ( methylation_pos[context] / ( methylation_neg[context] + methylation_pos[context]) ).T.to_pickle(f'{output_alias}.{context}.beta.pickle.gz')
