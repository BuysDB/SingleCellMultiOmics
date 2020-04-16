#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pysam


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      description='Filter VCF file for the allele frequency of the reference allele. Only records with af>=args.af and af<=(1-args.af) are selected')

    argparser.add_argument('vcf', type=str)
    argparser.add_argument('-af', type=float, default=0.2, help='Minimum allele frequency of REF allele')
    argparser.add_argument('--V', action='store_true',  help='Invert selection')
    args = argparser.parse_args()
    try:
        with pysam.VariantFile(args.vcf) as v:
          print(v.header, end='')
          for record in v:
            for sample in record.samples:
                total = sum(record.samples[sample]['AD'])
                af = record.samples[sample]['AD'][0]/total
                if af>=args.af and af<=(1-args.af):
                    print(record, end='')
                elif args.V:
                    print(record, end='')
    except BrokenPipeError:
        pass
