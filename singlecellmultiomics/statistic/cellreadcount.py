#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections
import pandas as pd
import matplotlib
import numpy as np
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def readIsDuplicate(read):
    return (read.has_tag('RC') and read.get_tag('RC') > 1) or read.is_duplicate


class CellReadCount(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.read_counts = collections.Counter()
        self.molecule_counts = collections.Counter()

    def processRead(self, read):
        if not read.has_tag('SM'):
            return

        cell = read.get_tag('SM')

        if read.is_read2 and not read.is_proper_pair:
            return

        self.read_counts[cell] +=1

        if  not readIsDuplicate(read):
            return
        self.molecule_counts[cell] +=1

    def to_csv(self, path):
        pd.DataFrame({'reads':self.read_counts, 'umis':self.molecule_counts}).to_csv(path)

    def __repr__(self):
        return f'The average amount of reads is {np.mean(self.read_counts.values())}'

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots()
        ax.hist(self.read_counts.values())

        if title is not None:
            ax.set_title(title)

        ax.set_xlabel("# Reads")
        ax.set_ylabel("# Cells")
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()

        fig, ax = plt.subplots()
        ax.hist(self.molecule_counts.values())
        if title is not None:
            plt.title(title)

        ax.set_xlabel("# Molecules")
        ax.set_ylabel("# Cells")
        plt.tight_layout()
        plt.savefig(target_path.replace('.png', '.molecules.png'))
        plt.close()
