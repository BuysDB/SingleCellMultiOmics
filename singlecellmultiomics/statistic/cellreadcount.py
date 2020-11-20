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
import seaborn as sns


def readIsDuplicate(read):
    return (read.has_tag('RC') and read.get_tag('RC') > 1) or read.is_duplicate


class CellReadCount(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.read_counts = collections.Counter()
        self.molecule_counts = collections.Counter()

    def processRead(self, R1,R2=None):

        for read in [R1,R2]:
            if read is None:
                continue

            if not read.has_tag('SM'):
                continue

            cell = read.get_tag('SM')

            self.read_counts[cell] +=1

            if not read.is_duplicate:
                self.molecule_counts[cell] +=1
            break

    def to_csv(self, path):
        pd.DataFrame({'reads':self.read_counts, 'umis':self.molecule_counts}).to_csv(path)

    def __repr__(self):
        return f'The average amount of reads is {np.mean(list(self.read_counts.values()))}'

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots()
        print(self.read_counts)
        ax.hist(list(self.read_counts.values()), bins=25, zorder=1)

        if title is not None:
            ax.set_title(title)

        ax.set_xlabel("# Reads")
        ax.set_ylabel("# Cells")
        ax.grid(zorder=0)
        sns.despine()
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()

        fig, ax = plt.subplots()
        ax.hist(list(self.molecule_counts.values()), bins=25,zorder=1)
        ax.grid(zorder=0)
        sns.despine()
        if title is not None:
            plt.title(title)

        ax.set_xlabel("# Molecules")
        ax.set_ylabel("# Cells")
        plt.tight_layout()
        plt.savefig(target_path.replace('.png', '.molecules.png'))
        plt.close()
