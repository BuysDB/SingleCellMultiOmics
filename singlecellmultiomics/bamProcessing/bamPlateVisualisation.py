#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import singlecellmultiomics.pyutils as pyutils
from singlecellmultiomics.tagtools import tagtools
import pysamiterators.iterators as pysamIterators
import gzip
import pickle
import subprocess

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import string
import math
import singlecellmultiomics.modularDemultiplexer
TagDefinitions = singlecellmultiomics.modularDemultiplexer.TagDefinitions

argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Visualize single cell statistics on a plate plot')
argparser.add_argument('-o',  type=str, help="output plot folder path, every library will be visualised as a separate plate", default='./plots/')

argparser.add_argument('alignmentfiles',  type=str, nargs='*')
args = argparser.parse_args()

if not os.path.exists(args.o):
    os.makedirs(args.o)

# Visualize the following:
# PER LIBRARY / DEMUX method
# total fragments
# total fragments with correct site
# unique molecules

# 384 well format:

well2index = collections.defaultdict(dict)
index2well = collections.defaultdict(dict)
rows = string.ascii_uppercase[:16]
columns = list(range(1,25))
for ci in range(1,385):
    i = ci-1
    rowIndex = math.floor( i/len(columns) )
    row = rows[rowIndex]
    column = columns[i%len(columns)]
    well2index[384][(row,column)] = ci
    index2well[384][ci] = (row,column)



rawFragmentCount = collections.defaultdict( collections.Counter ) # (library, mux) -> cell -> counts
usableCount = collections.defaultdict( collections.Counter ) # (library, mux) -> cell -> counts
moleculeCount = collections.defaultdict( collections.Counter ) # (library, mux) -> cell -> counts

skipReasons = collections.Counter()
for alignmentFile in args.alignmentfiles:
    with pysam.AlignmentFile(alignmentFile) as f:
        for read in f:
            rawFragmentCount[(read.get_tag('LY'), read.get_tag('MX'))][read.get_tag('SM')] += 1

            if  read.get_tag('MX').startswith('CS2'):
                if read.has_tag('XT') or read.has_tag('EX'):
                    if read.is_read1: # We only count reads2
                        continue
                    usableCount[(read.get_tag('LY'), read.get_tag('MX'))][read.get_tag('SM')] += 1

                    if  read.has_tag('RC') and read.get_tag('RC')==1:
                        moleculeCount[(read.get_tag('LY'), read.get_tag('MX'))][read.get_tag('SM')] += 1
            else:

                if read.has_tag('DS'):
                    if not read.is_read1:
                        skipReasons['Not R1'] += 1
                        continue

                    usableCount[(read.get_tag('LY'), read.get_tag('MX'))][read.get_tag('SM')] += 1
                    if read.get_tag('RC')==1:
                        moleculeCount[(read.get_tag('LY'), read.get_tag('MX'))][read.get_tag('SM')] += 1
                else:
                    skipReasons['No DS'] += 1


# Visualize:
for data, name in [(rawFragmentCount,'raw_reads'),(usableCount,'usable_reads'),(moleculeCount, 'unique_molecules')]:

    for (library, mux), cellCounts in data.items():

        df = pd.DataFrame( {name:cellCounts} )
        if mux=='CS2C8U6':
            offset = 0
        else:
            offset = 1
        df['col'] =  [ index2well[384][(offset+int(x.rsplit('_')[-1]))][1]  for x in df.index]
        #df['row'] = [ index2well[384][(offset+int(x.rsplit('_')[-1]))][0]  for x in df.index]
        df['row'] = [ -rows.index(index2well[384][(offset+int(x.rsplit('_')[-1]))][0])  for x in df.index]

        df.plot.scatter(x='col',y='row',s=(df[name]/np.percentile(df[name],99)*200),
                       c=(0.2,0.2,0.5,0.9)
                       )
        plt.yticks(sorted(df['row'].unique())[::-1], sorted(rows) , rotation=0)
        plt.xticks(sorted(df['col'].unique()), sorted(columns) , rotation=0)
        plt.title('{name}, {mux}, library:{library')
        plt.savefig(args.o + f'/{name}_{mux}_{library}.png')

        plt.close()
