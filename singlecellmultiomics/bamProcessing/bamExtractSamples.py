#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pysam
import argparse
import sys
from singlecellmultiomics.bamProcessing import get_sample_to_read_group_dict
import collections
from singlecellmultiomics.utils.path import get_valid_filename

def extract_samples( input_bam_path, output_path, capture_samples, head=None ):
    """
    Extract samples from a bam file

    Args:
        input_bam_path(str) : path to bam file from which to extract data from specified samples

        output_path(str) : prefix of path to write bam files to

        capture_samples(dict) : dictionary of lists {group: [sample, sample, sample], }

        head(int) : write this amount of reads, then exit

    """
    bamFile = pysam.AlignmentFile(input_bam_path, "rb")
    header = bamFile.header.copy()

    # Write to multiple files:
    output_handles = {}
    print('Groups:')
    for group in capture_samples:

        output_handles[group] = pysam.AlignmentFile(output_path.replace('.bam',f'{group}.bam'), "wb", header=header)
        print(f'\t{group}\t{output_handles[group].filename.decode()}')

    sample2handle = {}
    for group,samples in capture_samples.items():
        for sample in samples:
            sample2handle[sample] = output_handles[group]

    written = 0
    for r in bamFile:

        try:
            sample2handle[r.get_tag('SM')].write(r)
        except KeyError:
            continue

        written += 1
        if head is not None and written >= head:
            break

    print(f'Filtering finished, {written} reads')
    for group,handle in output_handles.items():
        handle.close()
        try:
            os.system(f'samtools index {handle.filename.decode()}')
        except Exception as e:
            print(e)
            pass


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

    capture_samples = collections.defaultdict( set ) # group -> set( samples )
    with open(args.extract) as e:
        for line in e:
            parts = line.strip().split(None, 1)
            if len(parts)==0:
                continue

            if len(parts)==1:
                group = ''
                sample = parts[0]
            elif len(parts)==2:
                sample, group = parts
                group = get_valid_filename(group)
            else:
                print(f'Line with issue: {line} :: {parts}')
                raise ValueError("Please supply a file with one or two columns: [SAMPLE], or [SAMPLE]{tab}[GROUP]")
            capture_samples[group].add( parts[0] )

    extract_samples(args.bamfile, args.o, capture_samples, head=args.head )
