#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


class OversequencingHistogram(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

    def processRead(self, R1,R2=None):

        for read in [R1, R2]:
            if read is None:
                continue

            if read.has_tag('RC'):
                overseq = read.get_tag('RC')
                self.histogram[overseq] += 1
                if overseq > 1:
                    self.histogram[overseq - 1] -= 1

                # Compatibility with picard --TAG_DUPLICATE_SET_MEMBERS
                """
                java -jar `which picard.jar` MarkDuplicates I=sorted.bam O=marked_duplicates.bam M=marked_dup_metrics.txt TAG_DUPLICATE_SET_MEMBERS=1
                """
            elif read.has_tag('PG'):
                if read.has_tag('DS') and not read.is_duplicate:
                    self.histogram[read.get_tag('DS')] += 1
                else:
                    self.histogram[1] += 1

            break

    def __repr__(self):
        return f'The average oversequencing is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots()
        overseqRange = list(range(1, 20))

        ax.scatter(overseqRange, [self.histogram[x] for x in overseqRange])
        if title is not None:
            ax.set_title(title)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel("Amount of fragments associated with molecule")
        ax.set_ylabel("# Molecules")
        plt.tight_layout()
        plt.savefig(target_path)

        ax = plt.gca()
        ax.set_ylim(1, None)
        ax.set_yscale('log')
        plt.savefig(target_path.replace('.png', '.log.png'))

        plt.close()
