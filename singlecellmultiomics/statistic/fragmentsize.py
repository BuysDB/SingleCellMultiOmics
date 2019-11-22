#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def readIsDuplicate(read):
    return (read.has_tag('RC') and read.get_tag('RC') > 1) or read.is_duplicate


class FragmentSizeHistogram(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

        self.histogramReject = collections.Counter()
        self.histogramAccept = collections.Counter()

    def processRead(self, read):
        if not read.is_proper_pair or not read.is_paired or not read.is_read1 or readIsDuplicate(
                read):
            return

        mateStart = read.next_reference_start
        if mateStart is None:
            return
        readLen = read.infer_query_length()
        if readLen is None:
            return
        readA = (read.reference_start, read.reference_end)
        readB = (mateStart, mateStart + readLen)

        end = max(
            read.reference_start,
            read.reference_end,
            mateStart,
            mateStart + readLen)
        start = min(
            read.reference_start,
            read.reference_end,
            mateStart,
            mateStart + readLen)
        fragmentSize = end - start
        if fragmentSize > 1_000:
            return
        #print(fragmentSize, read.reference_start,  read.reference_end,mateStart,readLen  )
        self.histogram[fragmentSize] += 1

        if read.has_tag('DS') and not read.has_tag('RR'):
            self.histogramAccept[fragmentSize] += 1
        else:
            self.histogramReject[fragmentSize] += 1

    def __repr__(self):
        return f'The average fragment size is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots()
        ax.bar(
            list(
                self.histogram.keys()), list(
                self.histogram.values()), width=1)
        if title is not None:
            ax.set_title(title)

        ax.set_xlabel("Fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()

        fig, ax = plt.subplots()
        plt.bar(
            list(
                self.histogramReject.keys()), list(
                self.histogramReject.values()), width=1, color='r')
        if title is not None:
            plt.title(title)

        ax.set_xlabel("Rejected fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        plt.savefig(target_path.replace('.png', '.rejected.png'))
        plt.close()

        fig, ax = plt.subplots()
        ax.bar(
            list(
                self.histogramAccept.keys()), list(
                self.histogramAccept.values()), width=1, color='g')
        if title is not None:
            plt.title(title)

        ax.set_xlabel("Accepted fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        plt.savefig(target_path.replace('.png', '.accepted.png'))
        plt.close()
