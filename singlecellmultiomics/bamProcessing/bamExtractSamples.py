#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pysam
import argparse
import sys
from singlecellmultiomics.bamProcessing import get_sample_to_read_group_dict
import collections

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract samples/cells from bam file, writes to multiple bam files if multiple groups are supplied
    """)
    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('extract', help="""Sample file, a list of samples to extract. First column: sample name, second column (optional) group name""")
    argparser.add_argument(
        '-o',
        type=str,
        required=True,
        help='output bam path, to the end of the file name the group name is appended')
    argparser.add_argument(
        '--f',
        action='store_true',
        help='Force overwrite of existing files')
    argparser.add_argument('-head', type=int)
    args = argparser.parse_args()

    if os.path.exists(args.o) and not args.f:
        raise ValueError(
            f'The output file {args.o} already exists! Use --f or remove the file')

    assert args.o.endswith('.bam'), 'output file path should end with .bam'

    bamFile = pysam.AlignmentFile(args.bamfile, "rb")
    header = bamFile.header.copy()


    capture_samples = collections.defaultdict( set ) # group -> set( samples )
    with open(args.extract) as e:
        for line in e:
            parts = line.strip().split()

            if len(parts)==1:
                group = 0
                sample = parts[0]
            elif len(parts)==2:
                sample, group = parts
            else:
                raise ValueError("Please supply a file with one or two columns: [SAMPLE], or [SAMPLE]{tab}[GROUP]")
            capture_samples[group].add( parts[0] )

    # Write to multiple files:
    output_handles = {}
    print('Groups:')
    for group in capture_samples:

        output_handles[group] = pysam.AlignmentFile(args.o.replace('.bam',f'{group}.bam'), "wb", header=header)

        print(f'\t{group}\t{output_handles[group].filename.decode()}')

    written = 0
    for r in bamFile:
        if r.has_tag('SM'):
            sample = r.get_tag('SM')
            for group in capture_samples:
                if sample in capture_samples[group]:
                    # Write
                    output_handles[group].write(r)

            written += 1
            if args.head and written >= args.head:
                break

    print(f'Filtering finished, {written} reads')
    for group,handle in output_handles.items():
        handle.close()
        try:
            os.system(f'samtools index {handle.filename.decode()}')
        except Exception as e:
            print(e)
            pass
