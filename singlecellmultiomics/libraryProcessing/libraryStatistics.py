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




def select_bam_file(lookup):
    for l in lookup:
        if os.path.exists(l):
            return l
    return None

def select_fastq_file(lookup):
    for paths in lookup:
        if type(paths)==tuple:
            for path in paths:
                if os.path.exists(path):
                    return paths
        elif os.path.exists(paths):
            return (paths,)
    return None

if __name__=='__main__':
    argparser = argparse.ArgumentParser(
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     description='Obtain statistics from your libraries')
    argparser.add_argument('libraries',  type=str, nargs='*')
    argparser.add_argument('-head',  type=int)
    argparser.add_argument('-tagged_bam',  type=str, help='Alias of subpath to tagged bam file. For example /tagged/sorted.bam')
    argparser.add_argument('--v',  action='store_true')
    argparser.add_argument('--nort',  action='store_true')
    args = argparser.parse_args()

    for library in args.libraries:
        library_name = library
        if library.endswith('.bam'):
            # the library is a bam file..
            bamFile=os.path.abspath( library)
            library= os.path.dirname(os.path.abspath( bamFile) )
            library_name = os.path.basename( os.path.abspath( bamFile) )
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
            PlateStatistic(args),
            ScCHICLigation(args)

        ]

        demuxFastqFilesLookup = [
            (f'{library}/demultiplexedR1.fastq.gz', f'{library}/demultiplexedR2.fastq.gz'),
            (f'{library}/demultiplexedR1_val_1.fq.gz', f'{library}/demultiplexedR2_val_2.fq.gz'),
            (f'{library}/demultiplexedR1_val_1.fq', f'{library}/demultiplexedR2_val_2.fq')
        ]

        rejectFilesLookup = [
            (f'{library}/rejectsR1.fastq.gz', f'{library}/rejectsR2.fastq.gz'),
            (f'{library}/rejectsR1.fastq', f'{library}/rejectsR2.fastq'),
            (f'{library}/rejectsR1.fq.gz', f'{library}/rejectsR2.fq.gz'),
            (f'{library}/rejectsR1.fq', f'{library}/rejectsR2.fq'),
        ]

        taggedFilesLookup = [

            f'{library}/tagged/marked_duplicates.bam',
            f'{library}/tagged/resorted.featureCounts.bam',
            f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam',
            f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.bam',
            f'{library}/tagged/sorted.bam'
        ]
        if args.tagged_bam:
            taggedFilesLookup.append(library+'/'+args.tagged_bam)


        if 'cluster' in library:
            continue
        print(f'{Style.BRIGHT}Library {library}{Style.RESET_ALL}')
        # Check if the bam file is present
        if bamFile is None:
            bamFile = select_bam_file(taggedFilesLookup)


        statFile = f'{library}/statistics.pickle.gz'

        demuxFastqFiles = select_fastq_file(demuxFastqFilesLookup)
        rejectFastqFiles = select_fastq_file(rejectFilesLookup)

        print("Selected files:")
        print(demuxFastqFiles)
        print(rejectFastqFiles)
        print(bamFile)

        demuxReads = None
        rejectedReads=None
        if demuxFastqFiles is not None:
            firstMate = demuxFastqFiles[0]
            print(f'\tDemuxed > {firstMate}')
            if firstMate.endswith('.gz'):
                demuxReads = pyutils.wccountgz(firstMate)/4
            else:
                demuxReads = pyutils.wccount(firstMate)/4

        if rejectFastqFiles is not None:
            firstMate = rejectFastqFiles[0]
            print(f'\tRejects > {firstMate}')
            if firstMate.endswith('.gz'):
                rejectedReads = pyutils.wccountgz(firstMate)/4
            else:
                rejectedReads = pyutils.wccount(firstMate)/4

        if demuxReads is not None:
            rc.setRawDemuxCount(demuxReads, paired=True)

            if rejectedReads is not None:
                rc.setRawReadCount(rejectedReads+demuxReads, paired=True)

        if bamFile is not None and os.path.exists(bamFile):
            print(f'\tTagged > {bamFile}')
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
                if args.v:
                    print(e)

        with gzip.open(statFile,'wb') as f:
            pickle.dump(statDict, f)


        # Make plots:
        plot_dir = f'{library}/plots'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        for statistic in statistics:
            if not hasattr(statistic, 'plot'):
                print(f'Not making a plot for {statistic.__class__.__name__} as no plot method is defined')
                continue
            try:
                statistic.plot(f'{plot_dir}/{statistic.__class__.__name__}.png', title=library_name)
            except Exception as e:
                if args.v:
                    import traceback
                    traceback.print_exc()

        # Make tables:
        table_dir = f'{library}/tables'
        if not os.path.exists(table_dir):
            os.makedirs(table_dir)
        for statistic in statistics:
            if not hasattr(statistic, 'to_csv'):
                print(f'Not making a table for {statistic.__class__.__name__} as to_csv method is not defined')
                continue
            try:
                statistic.to_csv(f'{table_dir}/{statistic.__class__.__name__}_{library_name}.csv')
            except Exception as e:
                if args.v:
                    import traceback
                    traceback.print_exc()

        # Make RT reaction plot:
        if not args.nort:
            os.system(f"bamPlotRTstats.py {bamFile} -head 2_000_000 --notstrict -o {plot_dir}/RT_")
