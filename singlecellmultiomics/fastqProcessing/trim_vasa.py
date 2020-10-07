#!/usr/bin/env python
from more_itertools import chunked
import argparse
import gzip
import re

if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Trim vasa fastq files for homo-polymers
    """)


    argparser.add_argument('R2_fastq', metavar='R2_fastq', type=str)
    argparser.add_argument('R2_singletons_out', type=str)

    argparser.add_argument('-poly_length', type=int, default=10)
    argparser.add_argument('-min_read_len', type=int, default=20)
    args = argparser.parse_args()

    poly_A = args.poly_length*'A'
    poly_T = args.poly_length*'T'
    poly_G = args.poly_length*'G'

    poly_length = args.poly_length

    r2_trimmer = re.compile('[GA]*$')

    def trim_r2(header, sequence, comment, qualities, min_read_len ):

        start = sequence.find(poly_A)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        start = sequence.find(poly_G)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        # Trim any trailing A and G bases from the end and # Trim down 3 bases
        sequence = r2_trimmer.sub('',sequence)[:-3]
        qualities = qualities[:len(sequence)]


        return f'{header}{sequence.rstrip()}\n{comment}{qualities.rstrip()}\n', len(sequence)>=min_read_len

    def trim_r1(header, sequence, comment, qualities, min_read_len ):

        start = sequence.rfind(poly_T) # Trim poly T
        if start != -1:
            sequence, qualities = sequence[start+poly_length:], qualities[start+poly_length:]

        start = sequence.find(poly_G)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]


        return f'{header}{sequence.rstrip()}\n{comment}{qualities.rstrip()}\n', len(sequence)>=min_read_len

    with gzip.open(args.R2_fastq,'rt') as r2, \
         gzip.open(args.R2_singletons_out,'wt',compresslevel=1) as r2_single_out:

        for read2 in chunked(r2,4):
            (r2_o, valid_r2)  = trim_r2(*read2, args.min_read_len)
            if valid_r2:
                r2_single_out.write(r2_o)
