#!/hpc/hub_oudenaarden/mloenhout/miniconda3/envs/SCMO/bin/python3.6

# -*- coding: utf-8 -*-
import os
import pysam
import argparse
import sys

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Filter bam file for value of tag or multiple tags
    """)
    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('expression', help="""Filter expression. All reads which yield a boolean true value for this expression will be send to the output file. For example r.get_tag("SM")=="cell1"   or r.reference_start==10   . Take care that for ALL positions you should subtract 1 because counting will start at 0. See https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment for all attributes you can use for filtering. It is possible to filter for a compound expression. For example r.get_tag('BC')=='CAACTAGA' and r.reference_name=='1' """)
    argparser.add_argument(
        '-o',
        type=str,
        required=True,
        help='output bam path')
    argparser.add_argument(
        '--f',
        action='store_true',
        help='Force overwrite of existing files')
    argparser.add_argument('-head', type=int)
    args = argparser.parse_args()

    if os.path.exists(args.o) and not args.f:
        raise ValueError(
            f'The output file {args.o} already exists! Use --f or remove the file')

    bamFile = pysam.AlignmentFile(args.bamfile, "rb")
    header = bamFile.header.copy()
    outputFile = pysam.AlignmentFile(args.o, "wb", header=header)

    written = 0
    skipped = 0
    for r in bamFile:
        try:
            if eval(args.expression):
                outputFile.write(r)
                written += 1
        except KeyError:
            skipped +=1

        if args.head and written >= (args.head-skipped):
            break

    print(f'Filtering finished, {written} reads written, {skipped} reads skipped')
    outputFile.close()
    try:
        os.system(f'samtools index {args.o}')
    except Exception as e:
        pass
