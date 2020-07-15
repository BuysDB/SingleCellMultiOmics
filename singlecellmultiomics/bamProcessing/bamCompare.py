#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pysam
import argparse
import sys


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Compare molecule assignment in two bam files
    """)
    argparser.add_argument('origin', metavar='bamfile', type=str)
    argparser.add_argument('compare', metavar='bamfile', type=str)
    argparser.add_argument('contig', metavar='bamfile', type=str)
    argparser.add_argument('start', metavar='bamfile', type=int)
    argparser.add_argument('end', metavar='bamfile', type=int)

    args = argparser.parse_args()


    origin = pysam.AlignmentFile(args.origin)
    compare = pysam.AlignmentFile(args.compare)


    region = args.contig, args.start, args.end
    reads_origin = {}
    for read in origin.fetch(*region):
        reads_origin[(read.query_name, read.is_read1)] = read

    reads_compare = {}
    for read in compare.fetch(*region):
        reads_compare[(read.query_name, read.is_read1)] = read

    for (qname, is_read1), read in reads_origin.items():
        if not (qname, is_read1) in reads_compare:
            print('Read', (qname, is_read1), 'vanished' )
            print(read.get_tag('SM'), read.is_qcfail, read.is_duplicate)
            print(read)
