#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


class AlleleHistogram(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

    def processRead(self, R1,R2):

        for read in [R1,R2]:
            if read is None:
                continue
            if read.has_tag('DA'):
                self.histogram[read.get_tag('DA')] += 1

    def __repr__(self):
        rt = 'Allele observations:'
        for allele, obs in self.histogram.most_common():
            rt += f'{allele}\t:\t{obs}\n'
        return rt

    def __iter__(self):
        return iter(self.histogram.most_common())

    def plot(self, target_path, title=None):
        d = dict(self)
        fig, ax = plt.subplots()
        ax.scatter(list(d.keys()), list(d.values()))
        plt.subplots_adjust(hspace=1)
        ax.set_yscale('log')
        ax.set_ylabel('# Molecules')
        ax.set_xlabel('Times oversequenced')
        ax.set_xlim(0, 20.5)
        ax.set_ylim((1, None))

        if title is not None:
            plt.title(title)
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()
