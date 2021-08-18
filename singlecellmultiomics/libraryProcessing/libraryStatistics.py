#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from singlecellmultiomics.statistic import *
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import pysam
import collections
import argparse

from singlecellmultiomics.bamProcessing import bam_is_processed_by_program

from colorama import Fore
from colorama import Back
from colorama import Style
import singlecellmultiomics.pyutils as pyutils
from singlecellmultiomics.tagtools import tagtools
from pysamiterators import MatePairIteratorIncludingNonProper
import gzip
import pickle
import subprocess
from glob import glob

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def select_bam_file(lookup):
    for l in lookup:
        if os.path.exists(l):
            print(f'Found file at {l}')
            return l

    return None


def select_fastq_file(lookup):
    for paths in lookup:
        if isinstance(paths, tuple):
            for path in paths:
                if os.path.exists(path):
                    return paths
        elif os.path.exists(paths):
            return (paths,)
    return None


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Obtain statistics from your libraries')
    argparser.add_argument(
        'libraries',
        type=str,
        nargs='*',
        help="either a library structured folder, or the tagged BAM file")
    argparser.add_argument('-t', type=str, default='all-stats',
                           help="type of staistics to produce. options are \
            'demult-stats', 'meth-stats', 'chic-stats' and 'all-stats'")
    argparser.add_argument('-o', type=str, help="output file prefix")
    argparser.add_argument(
        '--plotsOnly',
        action='store_true',
        help="only make plots")
    argparser.add_argument(
        '--fatal',
        action='store_true',
        help="Fatal error on any issue")
    argparser.add_argument(
        '--sl',
        action='store_true',
        help="Show lookup paths")


    argparser.add_argument(
        '--tablesOnly',
        action='store_true',
        help="only make tables")

    argparser.add_argument('-head', type=int)
    argparser.add_argument(
        '-tagged_bam',
        type=str,
        help='Alias of subpath to tagged bam file. For example /tagged/sorted.bam')
    argparser.add_argument('--v', action='store_true')
    argparser.add_argument('--nort', action='store_true')
    argparser.add_argument('--nolorenz', action='store_true')
    args = argparser.parse_args()

    for library in args.libraries:
        if library.endswith('.bam'):
            # the library is a bam file..
            bamFile = os.path.abspath(library)
            library = os.path.dirname(os.path.abspath(bamFile))
            library_name = os.path.basename(os.path.abspath(bamFile))
            print("Bam file was supplied:")
            print(bamFile)
        else:
            print("A library was supplied, automatically detecting files ..")
            bamFile = None
            library_name = os.path.basename(library)
        rc = ReadCount(args)  # Is also mappability

        statistics = [
            rc,
            FragmentSizeHistogram(args),
            RejectionReasonHistogram(args),
            MappingQualityHistogram(args),
            OversequencingHistogram(args),
            CellReadCount(args)
        ]


        full_file_statistics = []

        if not args.nolorenz:
            full_file_statistics.append( Lorenz(args) )

        if(args.t in ['meth-stats', 'all-stats']):
            statistics.extend([
                MethylationContextHistogram(args),
                ConversionMatrix(args)
            ])

        if(args.t in ['chic-stats', 'all-stats']):
            statistics.extend([ScCHICLigation(args)])

        if(args.t in ['demult-stats', 'all-stats']):
            statistics.extend([
                TrimmingStats(args),
                AlleleHistogram(args),
                DataTypeHistogram(args),
                TagHistogram(args),
                PlateStatistic(args),
                PlateStatistic2(args)
            ])

        demuxFastqFilesLookup = [
            (f'{library}/demultiplexedR1.fastq.gz',
             f'{library}/demultiplexedR2.fastq.gz'),
            (f'{library}/demultiplexedR1_val_1.fq.gz',
             f'{library}/demultiplexedR2_val_2.fq.gz'),
            (f'{library}/demultiplexedR1_val_1.fq',
             f'{library}/demultiplexedR2_val_2.fq')]

        rejectFilesLookup = [
            (f'{library}/rejectsR1.fastq.gz', f'{library}/rejectsR2.fastq.gz'),
            (f'{library}/rejectsR1.fastq', f'{library}/rejectsR2.fastq'),
            (f'{library}/rejectsR1.fq.gz', f'{library}/rejectsR2.fq.gz'),
            (f'{library}/rejectsR1.fq', f'{library}/rejectsR2.fq'),
        ]

        taggedFilesLookup = [
            f'{library}/tagged.bam',
            f'{library}/tagged/tagged.bam',
            f'{library}/tagged/marked_duplicates.bam',
            f'{library}/tagged/resorted.featureCounts.bam',
            f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam',
            f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.bam',
            f'{library}/tagged/sorted.bam']
        if args.tagged_bam:
            taggedFilesLookup = [
                library + '/' + args.tagged_bam] + taggedFilesLookup

        if args.sl:
            print(f'{Style.BRIGHT}Demux file lookup paths{Style.RESET_ALL}')
            print(demuxFastqFilesLookup)
            print(f'{Style.BRIGHT}Reject fastq file lookup paths{Style.RESET_ALL}')
            print(rejectFilesLookup)
            print(f'{Style.BRIGHT}Tagged bam lookup paths{Style.RESET_ALL}')
            print(taggedFilesLookup)

        if 'cluster' in library:
            continue
        print(f'{Style.BRIGHT}Library {library}{Style.RESET_ALL}')
        # Check if the bam file is present
        if bamFile is None:
            bamFile = select_bam_file(taggedFilesLookup)

        if bamFile is None:
            # Perform glob expansion
            bams = list(glob(f'{library}/*.bam'))+list(glob(f'{library}/*/*.bam'))
            for bam_path in bams:
                print(f"Trying {bam_path}",end="\t")
                try:
                    with pysam.AlignmentFile(bam_path) as a:
                        if bam_is_processed_by_program(a, program='bamtagmultiome'):
                            bamFile = bam_path
                            print("[TAGGED]")
                            break
                        else:
                            print("[NOT TAGGED]")
                except Exception as e:
                    print(f"[ERROR] {e}")

        if bamFile is None:
            print(f'{Fore.RED}BAM FILE MISSING {library}{Style.RESET_ALL}')
            exit()
        else:
            print(f'{Fore.GREEN}Bam file at {bamFile}{Style.RESET_ALL}')

        demuxFastqFiles = select_fastq_file(demuxFastqFilesLookup)
        rejectFastqFiles = select_fastq_file(rejectFilesLookup)

        print("Selected files:")
        print(f'demultiplexed reads: {demuxFastqFiles}')
        print(f'rejected reads: {rejectFastqFiles}')
        print(f'tagged bam: {bamFile}')

        demuxReads = None
        rejectedReads = None
        if demuxFastqFiles is not None:
            firstMate = demuxFastqFiles[0]
            print(f'\tDemuxed > {firstMate}')
            if firstMate.endswith('.gz'):
                demuxReads = pyutils.wccountgz(firstMate) / 4
            else:
                demuxReads = pyutils.wccount(firstMate) / 4

        if rejectFastqFiles is not None:
            firstMate = rejectFastqFiles[0]
            print(f'\tRejects > {firstMate}')
            if firstMate.endswith('.gz'):
                rejectedReads = pyutils.wccountgz(firstMate) / 4
            else:
                rejectedReads = pyutils.wccount(firstMate) / 4

        if demuxReads is not None:
            rc.setRawDemuxCount(demuxReads, paired=True)

            if rejectedReads is not None:
                rc.setRawReadCount(rejectedReads + demuxReads, paired=True)

        if bamFile is not None and os.path.exists(bamFile):
            print(f'\tTagged > {bamFile}')
            with pysam.AlignmentFile(bamFile) as f:

                for i, (R1,R2) in enumerate(MatePairIteratorIncludingNonProper(f)):
                    for statistic in statistics:
                        statistic.processRead(R1,R2)
                    if args.head is not None and i >= (args.head - 1):
                        break
        else:
            print(f'Did not find a bam file at {bamFile}')

        statDict = {}

        if os.path.exists(
                f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam'):
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

        for statistic in statistics:
            try:
                print(f'\t{statistic.__class__.__name__}')
                print(f'\t\t{statistic}\n')
                statDict[statistic.__class__.__name__] = dict(statistic)
                print(dict(statistic))
            except Exception as e:
                if args.v:
                    print(e)
                if args.fatal:
                    raise

        if bamFile is not None:
            for statistic in full_file_statistics:
                try:
                    with pysam.AlignmentFile(bamFile ) as alignments:
                        statistic.process_file(alignments)
                except Exception as e:
                    if args.v:
                        print(e)
                    if args.fatal:
                        raise



        # Make plots:
        if args.o is None:
            plot_dir = f'{library}/plots'
            table_dir = f'{library}/tables'
            statFile = f'{library}/statistics.pickle.gz'
        else:
            plot_dir = f'{args.o}/plots'
            table_dir = f'{args.o}/tables'
            statFile = f'{args.o}/statistics.pickle.gz'

        if not args.tablesOnly:

            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
            for statistic in statistics + full_file_statistics:
                if not hasattr(statistic, 'plot'):
                    print(
                        f'Not making a plot for {statistic.__class__.__name__} as no plot method is defined')
                    continue
                try:
                    statistic.plot(
                        f'{plot_dir}/{statistic.__class__.__name__}.png',
                        title=library_name)
                except Exception as e:
                    if args.v:
                        import traceback
                        traceback.print_exc()

        # Make tables:
        if not args.plotsOnly:
            with gzip.open(statFile, 'wb') as f:
                pickle.dump(statDict, f)
            if os.path.exists(statFile):
                with gzip.open(statFile, 'rb') as f:
                    try:
                        statDict.update(pickle.load(f))
                    except Exception as e:
                        if args.fatal:
                            raise


            if not os.path.exists(table_dir):
                os.makedirs(table_dir)
            for statistic in statistics + full_file_statistics:
                if not hasattr(statistic, 'to_csv'):
                    print(
                        f'Not making a table for {statistic.__class__.__name__} as to_csv method is not defined')
                    continue
                try:
                    statistic.to_csv(
                        f'{table_dir}/{statistic.__class__.__name__}_{library_name}.csv')
                except Exception as e:
                    if args.fatal:
                        raise
                    if args.v:
                        import traceback
                        traceback.print_exc()


        # Make RT reaction plot:
        if bamFile is not None and os.path.exists(bamFile):
            if not args.nort:
                os.system(
                    f"bamPlotRTstats.py {bamFile} -head 2_000_000 --notstrict -o {plot_dir}/RT_")
