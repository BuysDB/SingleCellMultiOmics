#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from singlecellmultiomics.utils import create_fasta_dict_file


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Create fasta dict file')

    argparser.add_argument('fastafile', type=str)

    args = argparser.parse_args()

    dict_path = create_fasta_dict_file(args.fastafile)
