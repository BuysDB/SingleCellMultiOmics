#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse

from colorama import Fore
from colorama import Back
from colorama import Style
import singlecellmultiomics.pyutils as pyutils
from singlecellmultiomics.tagtools import tagtools
import pysamiterators.iterators as pysamIterators
import gzip
import pickle
import subprocess

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
from singlecellmultiomics.statistic import *


argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Obtain statistics from your libraries')
argparser.add_argument('libraries',  type=str, nargs='*')
argparser.add_argument('-head',  type=int)
args = argparser.parse_args()

for library in args.libraries:
    if library.endswith('.bam'):
        # the library is a bam file..
        bamFile=os.path.abspath( library)
        library= os.path.dirname(os.path.abspath( bamFile) )
        print("Bam file was supplied:")
        print(bamFile)
    else:
        bamFile=None
    rc = ReadCount(args) # Is also mappability

    statistics = [
        rc,
        MethylationContextHistogram(args),
        MappingQualityHistogram(args),
        OversequencingHistogram(args),
        FragmentSizeHistogram(args),
        TrimmingStats(args),
        AlleleHistogram(args),
        RejectionReasonHistogram(args),
        DataTypeHistogram(args),
        TagHistogram(args),
        PlateStatistic(args)
    ]

    if 'cluster' in library:
        continue
    print(f'{Style.BRIGHT}Library {library}{Style.RESET_ALL}')
    # Check if the bam file is present
    if bamFile is None:
        if os.path.exists(f'{library}/tagged/sorted.bam'):
            bamFile = f'{library}/tagged/sorted.bam'
        if os.path.exists(f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam'):
            bamFile = f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam'
        if os.path.exists(f'{library}/tagged/resorted.featureCounts.bam'):
            bamFile = f'{library}/tagged/resorted.featureCounts.bam'

    statFile = f'{library}/statistics.pickle.gz'

    demuxFastqFiles = (f'{library}/demultiplexedR1.fastq.gz', f'{library}/demultiplexedR2.fastq.gz')
    demuxFastqFiles_alt = (f'{library}/demultiplexedR1_val_1.fq', f'{library}/demultiplexedR2_val_2.fq')
    rejectFastqFiles =  (f'{library}/rejectsR1.fastq.gz', f'{library}/rejectsR2.fastq.gz')


    rejectedReads =  pyutils.wccountgz(rejectFastqFiles[0])/4


    if os.path.exists(demuxFastqFiles[0]) and os.path.exists(demuxFastqFiles[1]):
        demuxReads = pyutils.wccountgz(demuxFastqFiles[0])/4
        # Perform fastq line count
        print(f'\t> {demuxFastqFiles[0]}')
        rc.setRawReadCount(rejectedReads+demuxReads, paired=True)

    else:
        if os.path.exists(demuxFastqFiles_alt[0]):
            demuxReads = pyutils.wccountgz(demuxFastqFiles_alt[0])/4
            rc.setRawReadCount(rejectedReads+demuxReads, paired=True)
        print(f'Did not find demultiplexed fastq files at {demuxFastqFiles[0]} {demuxFastqFiles[1]}' )


    if os.path.exists(demuxFastqFiles[0]) and os.path.exists(demuxFastqFiles[1]):
        # Perform fastq line count
        print(f'\t> {demuxFastqFiles[0]}')
        #reads = pyutils.wccountgz(demuxFastqFiles[0])/4
        rc.setRawDemuxCount(demuxReads, paired=True)
    else:
        print(f'Did not find demultiplexed fastq files at {demuxFastqFiles[0]} {demuxFastqFiles[1]}' )

    if bamFile is not None and os.path.exists(bamFile):
        print(f'\t> {bamFile}')
        with pysam.AlignmentFile(bamFile) as f:
            for i,read in enumerate(f):
                for statistic in statistics:
                    statistic.processRead(read)
                if args.head is not None and i>=(args.head-1):
                    break
    else:
        print(f'Did not find a bam file at {bamFile}' )

    statDict = {}

    if os.path.exists(f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam'):
        cmd = f'samtools view {bamFile} -F 4 -f 64 | cut -f 1 | sort | uniq | wc -l'
        out = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True
                             ).communicate()[0]
        read1mapped = int(out.partition(b' ')[0])
        cmd = f'samtools view {bamFile} -F 4 -f 128 | cut -f 1 | sort | uniq | wc -l'
        out = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True
                             ).communicate()[0]
        read2mapped = int(out.partition(b' ')[0])
        print(read1mapped,read2mapped)
        rc.totalMappedReads['R1'] = read1mapped
        rc.totalMappedReads['R2'] = read2mapped
        # Deduplicated reads have RC:i:1 set, -f 64 selects for R2
        cmd = f'samtools view {bamFile} -F 4 -f 64 | grep RC:i:1 | cut -f 1 | sort | uniq | wc -l'
        out = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True
                             ).communicate()[0]
        read1mappeddedup = int(out.partition(b' ')[0])
        # Deduplicated reads have RC:i:1 set, -f 128 selects for R2
        cmd = f'samtools view {bamFile} -F 4 -f 128 | grep RC:i:1 |  cut -f 1 | sort | uniq | wc -l'
        out = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True
                             ).communicate()[0]
        read2mappeddedup = int(out.partition(b' ')[0])

        rc.totalDedupReads['R1'] = read1mappeddedup
        rc.totalDedupReads['R2'] = read2mappeddedup




    if os.path.exists(statFile):
        with gzip.open(statFile,'rb') as f:
            try:
                statDict.update(pickle.load(f))
            except Exception as e:
                pass


    for statistic in statistics:
        try:
            print(f'\t{statistic.__class__.__name__}')
            print(f'\t\t{statistic}\n')
            statDict[statistic.__class__.__name__] = dict(statistic)
            print(dict(statistic))
        except Exception as e:
            print(e)

    with gzip.open(statFile,'wb') as f:
        pickle.dump(statDict, f)


    # Make plots:
    plot_dir = f'{library}/plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    for statistic in statistics:
        try:
            statistic.plot(f'{plot_dir}/{statistic.__class__.__name__}.png', title=library)
        except Exception as e:
            import traceback
            traceback.print_exc()
    # Make RT reaction plot:
    os.system(f"bamPlotRTstats.py {bamFile} -head 2_000_000 --notstrict -o {plot_dir}/RT_")
