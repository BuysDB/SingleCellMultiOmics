#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .statistic import Statistic
import singlecellmultiomics.pyutils as pyutils
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import collections

import pandas as pd
import numpy as np


class ReadCount(Statistic):
    def __init__(self,args):
        Statistic.__init__(self, args)
        self.totalMappedReads = collections.Counter()
        self.rawReadCount = None
        self.unmappedReads = collections.Counter()
        self.totalDedupReads = collections.Counter()
        self.totalAssignedSiteReads = collections.Counter()
        self.rejectionReasons = collections.Counter()
        self.demuxReadCount=0
        self.rawReadCount=0

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
