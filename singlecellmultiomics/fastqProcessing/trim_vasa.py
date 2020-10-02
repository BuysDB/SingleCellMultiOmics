#!/usr/bin/env python
from more_itertools import chunked
import argparse
import gzip

if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract methylation calls from bam file
    """)
    argparser.add_argument('R1_fastq', metavar='R1_fastq', type=str)
    argparser.add_argument('R2_fastq', metavar='R2_fastq', type=str)
    argparser.add_argument('R1_fastq_out', type=str)
    argparser.add_argument('R2_fastq_out', type=str)
    argparser.add_argument('R1_singletons_out', type=str)
    argparser.add_argument('-poly_length', type=int, default=10)
    args = argparser.parse_args()

    poly_A = args.poly_length*'A'
    poly_T = args.poly_length*'T'
    poly_G = args.poly_length*'G'

    poly_length = args.poly_length

    def trim_r2(header, sequence, comment, qualities ):

        start = sequence.find(poly_A)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        start = sequence.find(poly_G)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        return f'{header}{sequence.rstrip()}\n{comment}{qualities.rstrip()}\n', len(sequence)>=5

    def trim_r1(header, sequence, comment, qualities ):

        start = sequence.rfind(poly_T) # Trim poly T
        if start != -1:
            sequence, qualities = sequence[start+poly_length:], qualities[start+poly_length:]

        start = sequence.find(poly_G)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]


        return f'{header}{sequence.rstrip()}\n{comment}{qualities.rstrip()}\n', len(sequence)>=5


    with gzip.open(args.R1_fastq,'rt') as r1, \
         gzip.open(args.R2_fastq,'rt') as r2, \
         gzip.open(args.R1_fastq_out,'wt',compresslevel=1) as r1_out, \
         gzip.open(args.R1_singletons_out,'wt',compresslevel=1) as r1_single_out, \
         gzip.open(args.R2_fastq_out,'wt',compresslevel=1) as r2_out :
        for read1, read2 in zip(chunked(r1,4), chunked(r2,4)):
            (r1_o, valid_r1), (r2_o, valid_r2) = trim_r1(*read1), trim_r2(*read2)
            if valid_r1 and valid_r2:
                r1_out.write(r1_o), r2_out.write(r2_o)
            elif valid_r1:
                r1_single_out.write(r1_o)
