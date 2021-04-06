#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pyBigWig
from multiprocessing import Pool
from glob import glob
import os
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from singlecellmultiomics.methylation import sort_methylation_tabfile
from singlecellmultiomics.utils import create_fasta_dict_file


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Convert methylation calls from TAPS-tabulator to a methylation (big)bed file')

    argparser.add_argument('tabulated_file', type=str)
    argparser.add_argument(
        '-ref',
        type=str,
        help='path to reference fasta file, only required for bigbed files')

    argparser.add_argument(
        '-allele',
        type=str,
        help='Only generate calls for the specified allele')

    argparser.add_argument(
        '-o',
        type=str,
        help='Output path, ends in .bed to get a bed file and .bb to get a compressed bigbed file')

    args = argparser.parse_args()

    tabulated_file = args.tabulated_file
    # First we need to sort the callfile
    if tabulated_file.endswith('.sorted.gz'):
        pass
    else:
        print("Input seems not to be sorted. Performing sort")
        pathout = path += '.sorted.gz'
        sort_methylation_tabfile(path, pathout, True)
        tabulated_file = pathout

    if args.o.endswith('.bb'):
        print('exporting to big-bed file')
        if args.ref is None:
            raise ValueError('-ref is required for bigbed output')

        dict_path = create_fasta_dict_file(args.ref)
        # We need to write the bed file temporarily
        output_path = args.o.replace('.bb','.bed')
        methylation_tabfile_to_bed(tabulated_file, output_path)

        # Perform bigbed conversion:
        os.system(f'bedToBigBed {output_path} {dict_path} {args.o} -type=bed8+3')
        if os.path.exists(args.o):
            os.remove(output_path)
            print('Completed')
        else:
            print('bedToBigBed failed, is bedToBigBed installed ?')
    else:
        print('exporting to bed file')
        output_path = args.o
        methylation_tabfile_to_bed(tabulated_file, args.o)
