#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pysam
import argparse
import sys
import collections
from multiprocessing import Pool
from singlecellmultiomics.utils.path import get_valid_filename

def index_bam(path):
    try:
        pysam.index(path)
    except Exception as e:

        print(e)

def split_bam_by_tag( input_bam_path, output_prefix, tag, head=None, max_handles=10, skip=None ):
    """
    Split bam file by a tag value

    Args:
        input_bam_path(str) : path to bam file from which to extract data from specified samples

        output_prefix(str) : prefix of path to write bam files to

        tag(str) : tag to split file on

        head(int) : write this amount of reads, then exit

    """
    bamFile = pysam.AlignmentFile(input_bam_path, "rb")
    header = bamFile.header.copy()

    # Write to multiple files:
    output_handles = {}

    done = set()
    waiting = set()

    print('Started writing', end='')
    written = 0
    for r in bamFile:
        if not r.has_tag(tag):
            continue
        value = str(r.get_tag(tag))

        # Clean up value to get proper file name:
        value = get_valid_filename(value)

        if value in skip or value in waiting:
            continue


        if not value in output_handles:
            if len(output_handles)>=max_handles:
                waiting.add(value)
                continue
            print(f'\rOpened bam file for {value}            ', end='')
            output_handles[value] = pysam.AlignmentFile(f'{output_prefix}.{value}.bam', "wb", header=header)

        output_handles[value].write(r)
        written+=1
        if written%100_000==0:
            print(f'\rWrote {written} reads                      ',end='')
        if head is not None and written>=head:
            break


    to_index = []
    for value,handle in output_handles.items():
        done.add(value)
        handle.close()
        to_index.append(handle.filename.decode())


    # Index using multiple threads:
    print(f'\nIndexing {len(to_index)} bam files')
    with Pool(10) as workers:
        for result in workers.imap_unordered(index_bam, to_index):
            pass

    return done, waiting

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Split bam file by a tag, writes reads with a different value for the specified tag to multiple bam files
    """)

    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('-max_handles', default=400, type=int)
    argparser.add_argument('tag', help="tag to use for splitting")
    argparser.add_argument(
        '-o',
        type=str,
        required=True,
        help='output bam prefix, to the end of the file name the tag value is appended')

    argparser.add_argument('-head', type=int)
    args = argparser.parse_args()

    skip= set()

    waiting = set([0]) #
    iteration = 1
    while len(waiting)>0:
        print(f'Iteration {iteration}')
        done, waiting = split_bam_by_tag(args.bamfile, args.o, args.tag, head=args.head, max_handles=args.max_handles,skip=skip )
        print('Wrote bam files for tag values:')
        for d in done:
            print(f'\t{d}')
        skip.update(done)
        iteration+=1
    print(f'All done, wrote {len(skip)} bam files')
