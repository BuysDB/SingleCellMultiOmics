#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import pandas as pd
import singlecellmultiomics
import singlecellmultiomics.modularDemultiplexer
TagDefinitions = singlecellmultiomics.modularDemultiplexer.TagDefinitions


argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Obtain duplication rate on the fly')

argparser.add_argument('alignmentfile',  type=str)

args = argparser.parse_args()


molecule_count = 0
read_count = 0
f = pysam.AlignmentFile(args.alignmentfile)
for i, molecule in enumerate(singlecellmultiomics.molecule.MoleculeIterator(
    f,
    fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
    contig='chr1'
    )):

    molecule_count+=1
    read_count += len(molecule)

    if i%10_000==0:
        print(f'\r{read_count/molecule_count}             ',end='')
