#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


class ScCHICLigation():
    def __init__(self, args):
        # cell -> { A_start: count, total_cuts: count }
        self.per_cell_a_obs = collections.defaultdict(collections.Counter)
        # cell -> { TA_start: count, total_cuts: count }
        self.per_cell_ta_obs = collections.defaultdict(collections.Counter)

    def processRead(self, R1,R2):

        if R1 is None:
            return
        read = R1
        if read.has_tag('RZ') and not read.is_duplicate and read.is_read1:
            sample = read.get_tag('SM')
            first = read.get_tag('RZ')[0]
            if read.get_tag('RZ') == 'TA':
                self.per_cell_ta_obs[sample]['TA_start'] += 1
            if first == 'A':
                self.per_cell_a_obs[sample]['A_start'] += 1
            self.per_cell_ta_obs[sample]['total'] += 1
            self.per_cell_a_obs[sample]['total'] += 1

    def __repr__(self):
        return 'ScCHICLigation: no description'



    def __iter__(self):
        for cell, cell_data in self.per_cell_ta_obs.items():
            yield cell_data['total'],  cell_data['TA_start'] / cell_data['total']

    def plot(self, target_path, title=None):

        ########### TA ###########
        fig, ax = plt.subplots(figsize=(4, 4))

        x = []
        y = []
        for cell, cell_data in self.per_cell_ta_obs.items():
            x.append(cell_data['total'])
            y.append(cell_data['TA_start'] / cell_data['total'])

        ax.scatter(x, y, s=3,c='k')
        ax.set_xscale('log')
        if title is not None:
            ax.set_title(title)

        ax.set_ylabel("Fraction unique cuts starting with TA")
        ax.set_xlabel("# Molecules")
        ax.set_xlim(1, None)
        ax.set_ylim(-0.1, 1.05)
        sns.despine()
        plt.tight_layout()
        plt.savefig(target_path.replace('.png', '.TA.png'))
        plt.close()

        ########### A ###########
        fig, ax = plt.subplots(figsize=(4, 4))

        x = []
        y = []
        for cell, cell_data in self.per_cell_ta_obs.items():
            x.append(cell_data['total'])
            y.append(cell_data['A_start'] / cell_data['total'])

        ax.scatter(x, y, s=3,c='k')
        ax.set_xscale('log')
        if title is not None:
            ax.set_title(title)

        ax.set_ylabel("Fraction unique cuts starting with A")
        ax.set_xlabel("# Molecules")
        ax.set_xlim(1, None)
        ax.set_ylim(-0.1, 1.05)
        plt.tight_layout()
        sns.despine()
        plt.savefig(target_path.replace('.png', '.A.png'))
        plt.close()

    def to_csv(self, path):
        pd.DataFrame(
            self.per_cell_ta_obs).sort_index().to_csv(
            path.replace(
                '.csv',
                'TA_obs_per_cell.csv'))
