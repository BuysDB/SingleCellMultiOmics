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
import matplotlib.lines as mlines
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import string
import math
import singlecellmultiomics.modularDemultiplexer
TagDefinitions = singlecellmultiomics.modularDemultiplexer.TagDefinitions
import matplotlib.patheffects as path_effects

def human_readable( value, targetDigits=2,fp=0):

    #Float:
    if value<1 and value>0:
        return('%.2f' % value )

    if value == 0.0:
        return('0')

    baseId = int(math.floor( math.log10(float(value))/3.0 ))
    suffix = ""
    if baseId==0:
        sVal =  str(round(value,targetDigits))
        if len(sVal)>targetDigits and sVal.find('.'):
            sVal = sVal.split('.')[0]

    elif baseId>0:

        sStrD = max(0,targetDigits-len(str( '{:.0f}'.format((value/(math.pow(10,baseId*3)))) )))


        sVal = ('{:.%sf}' % min(fp, sStrD)).format((value/(math.pow(10,baseId*3))))
        suffix = 'kMGTYZ'[baseId-1]
    else:

        sStrD = max(0,targetDigits-len(str( '{:.0f}'.format((value*(math.pow(10,-baseId*3)))) )))
        sVal = ('{:.%sf}' %  min(fp, sStrD)).format((value*(math.pow(10,-baseId*3))))
        suffix = 'mnpf'[-baseId-1]

        if len(sVal)+1>targetDigits:
            # :(
            sVal = str(round(value,fp))[1:]
            suffix = ''


    return('%s%s' % (sVal,suffix))

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
        df['size'] = (df[name]/np.percentile(df[name],99)*200)
        df.plot.scatter(x='col',y='row',s=df['size'],
                       c=[(0.2,0.2,0.5,0.9)]
                       )

        ax = plt.gca()
        for ii,row in df.iterrows():
            if row[name]>0 and (row[name]<np.percentile(df[name],5) or row[name]>np.percentile(df[name],95)):
                text = ax.annotate( human_readable(int(row[name])), ( row['col'], row['row']),
                ha='center',va='center_baseline',color='w',size=7)
                text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                       path_effects.Normal()])

        plt.yticks(sorted(df['row'].unique())[::-1], sorted(rows) , rotation=0)
        plt.xticks(sorted(df['col'].unique()), sorted(columns) , rotation=0)
        plt.title(fr'{name} with ${mux}$ adapter'+f'\n{library}')

        # Create legend:
        #ld = []
        #for x in np.linspace(1, max(df[name]), 4):
    #        size = (x/np.percentile(df[name],99))*200
        #    ld.append( mlines.Line2D([], [], color='blue', marker='.', linestyle='None',
        ##                  markersize=np.sqrt(size), label=f'{int(x)}:{size}'))
        #plt.legend(handles=ld)

        plt.savefig(args.o + f'/{name}_{mux}_{library}.png')

        plt.close()
