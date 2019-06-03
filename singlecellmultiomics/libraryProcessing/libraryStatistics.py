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



argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Obtain statistics from your libraries')
argparser.add_argument('libraries',  type=str, nargs='*')
argparser.add_argument('-head',  type=int)
args = argparser.parse_args()


def readIsDuplicate(read):
    return read.has_tag('RC') and read.get_tag('RC')>1

class Statistic(object):

    def __init__(self, args):
        self.args = args

    def processRead(self,read):
        pass

    def __repr__(self):
        return 'dummy'


class ReadCount(Statistic):
    def __init__(self,args):
        Statistic.__init__(self, args)
        self.totalMappedReads = collections.Counter()
        self.rawReadCount = None
        self.unmappedReads = collections.Counter()
        self.totalDedupReads = collections.Counter()
        self.totalAssignedSiteReads = collections.Counter()
        self.rejectionReasons = collections.Counter()

    def plot(self, target_path, title=None):
        df = pd.DataFrame.from_dict(dict(self))
        df['Raw reads']/=2 # They are paired and only present one time in the dict but are expanded by pandas
        df['Demultiplexed reads']/=2 # same
        df = df[['Raw reads', 'Demultiplexed reads', 'Mapped reads','AssignedSiteReads','Deduplicated reads']] #,'UnmappedReads']]
        df.plot.bar(figsize=(10,4)).legend(bbox_to_anchor=(1, 0.98))
        if title is not None:
            plt.title(title)

        plt.subplots_adjust(right=0.6)
        plt.savefig(target_path)

        ax = plt.gca()
        ax.set_yscale('log')
        plt.savefig(target_path.replace('.png','.log.png'))

        plt.close()

    def processRead(self,read):

        # Count every read only once
        if read.is_supplementary or read.is_secondary:
            return


        if not read.is_unmapped:

            if read.is_read1:
                self.totalMappedReads['R1']+=1
            elif read.is_read2:
                self.totalMappedReads['R2']+=1
            else:
                self.totalMappedReads['R?']+=1
        else:
            if read.is_read1:
                self.unmappedReads['R1']+=1
            elif read.is_read2:
                self.unmappedReads['R2']+=1
            else:
                self.unmappedReads['R?']+=1

        if read.has_tag('DS'):

            if read.is_read1:
                self.totalAssignedSiteReads['R1']+=1
            else:
                self.totalAssignedSiteReads['R2']+=1

        if read.has_tag('RC') and read.get_tag('RC')==1:
            if  not read.is_secondary:
                if read.is_read1:
                    self.totalDedupReads['R1']+=1
                else:
                    self.totalDedupReads['R2']+=1


    def setRawReadCount(self, readCount, paired=True):
        self.rawReadCount = readCount *(2 if paired else 1)

    def setRawDemuxCount(self, readCount, paired=True):
        self.demuxReadCount = readCount *(2 if paired else 1)

    def mappability(self):
        if self.rawReadCount>0:
            return sum(self.totalMappedReads.values()) / self.rawReadCount
        return None
    def __repr__(self):
        return f"Input:{int(self.rawReadCount)} raw reads, demultiplexed:{self.demuxReadCount}, final mapped: {self.totalMappedReads['R1']} R1 reads, {self.totalMappedReads['R1']} R2 reads, {self.totalMappedReads['R?']} unknown pairing reads, mappability:{self.mappability() if self.mappability() else 'Cannot be determined'}"

    def __iter__(self):
        yield 'Raw reads', self.rawReadCount
        yield 'Demultiplexed reads', self.demuxReadCount
        yield 'Mapped reads', self.totalMappedReads
        yield 'UnmappedReads', self.unmappedReads
        yield 'AssignedSiteReads', self.totalAssignedSiteReads
        yield 'Deduplicated reads', self.totalDedupReads

class StatisticHistogram(Statistic):
    def __init__(self,args):
        Statistic.__init__(self, args)
        self.histogram = collections.Counter()

    def __repr__(self):
        return f'Mean {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(histogram)}'

    def __iter__(self):
        return iter(self.histogram.most_common())


class MappingQualityHistogram(StatisticHistogram):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()
    def processRead(self,read):
        self.histogram[read.mapping_quality]+=1
    def __repr__(self):
        return f'The average mapping quality is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'


class TrimmingStats(Statistic):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.totalFragmentsTrimmed = 0

    def processRead(self,read):
        if read.has_tag('a1') or read.has_tag('eB') or read.has_tag('A2') or read.has_tag('EB'):
            self.totalFragmentsTrimmed+=1

    def __repr__(self):
        return f'Trimmed fragments: {self.totalFragmentsTrimmed}'

    def __iter__(self):
        yield 'Trimmed fragments', self.totalFragmentsTrimmed

class OversequencingHistogram(StatisticHistogram):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()
    def processRead(self,read):
        if read.has_tag('RC'):
            overseq = read.get_tag('RC')
            self.histogram[overseq]+=1
            if overseq>1:
                self.histogram[overseq-1]-=1
    def __repr__(self):
        return f'The average oversequencing is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'


class AlleleHistogram(Statistic):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()
    def processRead(self,read):
        if read.has_tag('DA'):

            self.histogram[read.get_tag('DA')]+=1
    def __repr__(self):
        rt = 'Allele observations:'
        for allele, obs in self.histogram.most_common():
            rt += f'{allele}\t:\t{obs}\n'
        return rt
    def __iter__(self):
        return iter(self.histogram.most_common())

    def plot(self, target_path, title=None):
        d = dict(self)
        fig, ax  = plt.subplots()
        ax.scatter(list(d.keys()), list(d.values()))
        plt.subplots_adjust(hspace = 1)
        ax.set_yscale('log')
        ax.set_ylabel('# Molecules')
        ax.set_xlabel('Times oversequenced')
        ax.set_xlim(0,20.5)
        ax.set_ylim( (1,None))

        if title is not None:
            plt.title(title)
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()


class RejectionReasonHistogram(Statistic):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

    def processRead(self,read):
        if read.has_tag('RR'):
            self.histogram[read.get_tag('RR')]+=1
    def __repr__(self):
        rt = 'Rejection reasons:'
        for reason, obs in self.histogram.most_common():
            rt += f'{reason}\t:\t{obs}\n'
        return rt
    def __iter__(self):
        return iter(self.histogram.most_common())

    def plot(self, target_path, title=None):
        df = pd.DataFrame.from_dict({'Reason':dict(self)}).T

        df.plot.bar(figsize=(10,4)).legend(bbox_to_anchor=(1, 0.98))
        if title is not None:
            plt.title(title)

        plt.tight_layout()
        plt.subplots_adjust(right=0.6)
        plt.savefig(target_path)

        ax = plt.gca()
        ax.set_yscale('log')
        plt.savefig(target_path.replace('.png','.log.png'))
        plt.close()

class TagHistogram(Statistic):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

    def processRead(self,read):
        if read.has_tag('EX'):
            if read.get_tag('EX')=='Unassigned_NoFeatures':
                self.histogram["Not assigned to exon"]+=1
            elif read.get_tag('EX')=='Assigned':
                self.histogram["Assigned to exon"]+=1
            else:
                self.histogram["Unkown exon assignment"]+=1

        if read.has_tag('XS'):
            if read.get_tag('XS')=='Unassigned_NoFeatures':
                self.histogram["Not assigned to gene/intron"]+=1
            elif read.get_tag('XS')=='Assigned':
                self.histogram["Assigned to gene/intron"]+=1
            else:
                self.histogram["Unkown gene/intron assignment"]+=1

        if read.has_tag('Is'):
            self.histogram[f"Sequencer_{read.get_tag('Is')}"]+=1

        if read.has_tag('LY'):
            self.histogram[f"Library_{read.get_tag('LY')}"]+=1


    def __repr__(self):
        rt = 'Tag obs::'
        for reason, obs in self.histogram.most_common():
            rt += f'{reason}\t:\t{obs}\n'
        return rt
    def __iter__(self):
        return iter(self.histogram.most_common())

    def plot(self, target_path, title=None):
        df = pd.DataFrame.from_dict({'Tag':dict(self)}).T

        df.plot.bar(figsize=(10,4)).legend(bbox_to_anchor=(1, 0.98))
        if title is not None:
            plt.title(title)

        plt.tight_layout()
        plt.subplots_adjust(right=0.6)
        plt.savefig(target_path)

        ax = plt.gca()
        ax.set_yscale('log')
        plt.savefig(target_path.replace('.png','.log.png'))
        plt.close()

class DataTypeHistogram(Statistic):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

    def processRead(self,read):
        if read.has_tag('MX'):
            self.histogram[read.get_tag('MX')]+=1

    def __repr__(self):
        rt = 'Rejection reasons:'
        for dataType, obs in self.histogram.most_common():
            rt += f'{dataType}\t:\t{obs}\n'
        return rt
    def __iter__(self):
        return iter(self.histogram.most_common())

    def plot(self, target_path, title=None):
        df = pd.DataFrame.from_dict({'DataType':dict(self)}).T

        df.plot.bar(figsize=(10,4)).legend(bbox_to_anchor=(1, 0.98))
        if title is not None:
            plt.title(title)

        plt.tight_layout()
        plt.subplots_adjust(right=0.6)
        plt.savefig(target_path)

        ax = plt.gca()
        ax.set_yscale('log')
        plt.savefig(target_path.replace('.png','.log.png'))
        plt.close()



class FragmentSizeHistogram(StatisticHistogram):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

        self.histogramReject = collections.Counter()
        self.histogramAccept = collections.Counter()
    def processRead(self,read):
        if not read.is_paired or not read.is_read1 or readIsDuplicate(read):
            return

        mateStart = read.next_reference_start
        if mateStart is None :
            return
        readLen = read.infer_query_length()
        if readLen is None:
            return
        readA = (read.reference_start, read.reference_end)
        readB = (mateStart, mateStart+readLen)

        end = max(read.reference_start, read.reference_end, mateStart, mateStart+readLen)
        start = min(read.reference_start, read.reference_end, mateStart, mateStart+readLen)
        fragmentSize = end - start
        if fragmentSize>1_000:
            return
        #print(fragmentSize, read.reference_start,  read.reference_end,mateStart,readLen  )
        self.histogram[fragmentSize]+=1
        if read.has_tag('DS'):
            self.histogramAccept[fragmentSize]+=1
        else:
            self.histogramReject[fragmentSize]+=1

    def __repr__(self):
        return f'The average fragment size is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots()
        ax.bar( list(self.histogram.keys()), list(self.histogram.values() ),width=1)
        if title is not None:
            ax.set_title(title)

        ax.set_xlabel("Fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()

        fig, ax = plt.subplots()
        plt.bar( list(self.histogramReject.keys()), list(self.histogramReject.values() ),width=1,color ='r')
        if title is not None:
            plt.title(title)

        ax.set_xlabel("Rejected fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        plt.savefig(target_path.replace('.png','.rejected.png'))
        plt.close()


        fig, ax = plt.subplots()
        ax.bar( list(self.histogramAccept.keys()), list(self.histogramAccept.values() ),width=1,color ='g')
        if title is not None:
            plt.title(title)

        ax.set_xlabel("Accepted fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        plt.savefig(target_path.replace('.png','.accepted.png'))
        plt.close()


for library in args.libraries:
    bamFile=None
    rc = ReadCount(args) # Is also mappability
    statistics = [
        rc,
        MappingQualityHistogram(args),
        OversequencingHistogram(args),
        FragmentSizeHistogram(args),
        TrimmingStats(args),
        AlleleHistogram(args),
        RejectionReasonHistogram(args),
        DataTypeHistogram(args),
        TagHistogram(args)
    ]

    if 'cluster' in library:
        continue
    print(f'{Style.BRIGHT}Library {library}{Style.RESET_ALL}')
    # Check if the bam file is present
    if os.path.exists(f'{library}/tagged/sorted.bam'):
        bamFile = f'{library}/tagged/sorted.bam'
    if os.path.exists(f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam'):
        bamFile = f'{library}/tagged/STAR_mappedAligned.sortedByCoord.out.featureCounts.bam'
    if os.path.exists(f'{library}/tagged/resorted.featureCounts.bam'):
        bamFile = f'{library}/tagged/resorted.featureCounts.bam'
    statFile = f'{library}/statistics.pickle.gz'

    demuxFastqFiles = (f'{library}/demultiplexedR1.fastq.gz', f'{library}/demultiplexedR2.fastq.gz')
    rejectFastqFiles =  (f'{library}/rejectsR1.fastq.gz', f'{library}/rejectsR2.fastq.gz')

    demuxReads = pyutils.wccountgz(demuxFastqFiles[0])/4
    rejectedReads =  pyutils.wccountgz(rejectFastqFiles[0])/4

    if os.path.exists(demuxFastqFiles[0]) and os.path.exists(demuxFastqFiles[1]):
        # Perform fastq line count
        print(f'\t> {demuxFastqFiles[0]}')
        rc.setRawReadCount(rejectedReads+demuxReads, paired=True)

    else:
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
        print(f'\t{statistic.__class__.__name__}')
        print(f'\t\t{statistic}\n')
        statDict[statistic.__class__.__name__] = dict(statistic)
        print(dict(statistic))

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
