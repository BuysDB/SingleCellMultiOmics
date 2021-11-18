#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam
import sys
from datetime import datetime
import singlecellmultiomics
from singlecellmultiomics.fragment import Fragment
from singlecellmultiomics.bamProcessing.bamFunctions import get_read_group_from_read, sorted_bam_file, write_program_tag

def set_read_group( in_bam_path, out_bam_path, id:str, pl:str, lb:str, sm:str, pu:str, threads=4  ):
    """
    Set read group format of bam file

    Args:
        in_bam_path(str) : bam file to read

        in_bam_path(str) : bam file to write

        format(int) : formatting mode, see `singlecellmultiomics.fragment.Fragment.get_read_group`

        threads(int) : Amount of decompression threads for reading

    """

    """
    read_groups(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
            { 'ID': ID, 'LB':library,
            'PL':platform,
            'SM':sampleLib,
            'PU':readGroup }
    """

    read_groups =     {id:{ 'ID': id, 'LB':lb,
        'PL':pl,
        'SM':sm,
        'PU':pu }}

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
                rg_id = id
                read.set_tag('RG',rg_id)
                out.write(read)


if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Set read group sample platform unit etc.')
    argparser.add_argument('bamin', type=str)


    argparser.add_argument(
        '-id',
        type=str,
        required = True,
        help="Read group id")

    argparser.add_argument(
        '-pl',
        type=str,
        default='ILLUMINA',
        help="Read group platform")

    argparser.add_argument(
        '-pu',
        type=str,
        default='LANE',
        help="Read group platform unit")

    argparser.add_argument(
        '-sm',
        type=str,
        required=True,
        help="Read group sample name")

    argparser.add_argument(
        '-lb',
        type=str,
        default='mix',
        help="Library")


    argparser.add_argument(
        '-t',
        type=int,
        default=4,
        help="Threads")

    argparser.add_argument('-o', type=str, help="output bam file", required=True)
    args = argparser.parse_args()
    set_read_group(
        args.bamin,
        args.o,
        id = args.id,
        pl = args.pl,
        lb = args.lb,
        sm = args.sm,
        pu = args.pu,
        threads = args.t)
