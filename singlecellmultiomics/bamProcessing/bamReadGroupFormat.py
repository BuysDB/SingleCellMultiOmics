#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam
import sys
from datetime import datetime
import singlecellmultiomics
from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.bamProcessing.bamFunctions import get_read_group_from_read, sorted_bam_file, write_program_tag

def set_read_group_format( in_bam_path, out_bam_path, format, threads=4  ):
    """
    Set read group format of bam file

    Args:
        in_bam_path(str) : bam file to read

        in_bam_path(str) : bam file to write

        format(int) : formatting mode, see `singlecellmultiomics.fragment.Fragment.get_read_group`

        threads(int) : Amount of decompression threads for reading

    """
    read_groups = dict()  # Store unique read groups in this set
    """
    read_groups(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
            { 'ID': ID, 'LB':library,
            'PL':platform,
            'SM':sampleLib,
            'PU':readGroup }
    """
    with pysam.AlignmentFile(in_bam_path, threads = threads) as input_bam:

        input_header = input_bam.header.as_dict()

        # Write provenance information to BAM header
        write_program_tag(
            input_header,
            program_name='bamReadGroupFormat',
            command_line=" ".join(
                sys.argv),
            version=singlecellmultiomics.__version__,
            description=f'SingleCellMultiOmics read group formatting, executed at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

        with sorted_bam_file(out_bam_path, header=input_header, read_groups=read_groups, input_is_sorted=True) as out:
            print('Started writing')
            for read in input_bam:

                rg_id = get_read_group_from_read(read, format=format,  with_attr_dict=False  )
                if not rg_id in read_groups:
                    rg_id, rg_data = get_read_group_from_read(read, format=format, with_attr_dict=True  )
                    read_groups[rg_id] = rg_data

                read.set_tag('RG',rg_id)
                out.write(read)


if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Change format of read groups')
    argparser.add_argument('bamin', type=str)
    argparser.add_argument(
        '-format',
        type=int,
        default=0,
        help="Read group format code")

    argparser.add_argument(
        '-t',
        type=int,
        default=4,
        help="Threads")

    argparser.add_argument('-o', type=str, help="output bam file", required=True)
    args = argparser.parse_args()
    set_read_group_format(args.bamin, args.o, args.format, args.t)
