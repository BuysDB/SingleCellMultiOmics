#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import pandas as pd
import singlecellmultiomics
import singlecellmultiomics.molecule
import singlecellmultiomics.modularDemultiplexer
TagDefinitions = singlecellmultiomics.modularDemultiplexer.TagDefinitions


argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Obtain duplication rate on the fly')

argparser.add_argument('alignmentfile',  type=str)
argparser.add_argument('-u',  type=int, 10_000)

args = argparser.parse_args()


molecule_count = 0
read_count = 0
f = pysam.AlignmentFile(args.alignmentfile)

duprateDist = collections.Counter()



for i, molecule in enumerate(singlecellmultiomics.molecule.MoleculeIterator(
    f,
    fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment
    )):

    molecule_count+=1
    read_count += len(molecule)
    duprateDist[len(molecule)]+=1

    if i%args.u==0:
        print(f'\r{read_count/molecule_count}             ',end='')
