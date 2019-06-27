#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

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

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots()
        overseqRange = list(range(1,20))

        ax.scatter( overseqRange, [self.histogram[x] for x in overseqRange])
        if title is not None:
            ax.set_title(title)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel("Amount of fragments associated with molecule")
        ax.set_ylabel("# Molecules")
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()
