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

class ScCHICLigation():
    def __init__(self,args):
        self.per_cell_ta_obs = collections.defaultdict( collections.Counter ) # cell -> { A_start: count, total_cuts: count }

    def processRead(self,read):
        if read.has_tag('RZ') and not read.is_duplicate and read.is_read1:
            sample = read.get_tag('SM')
            first = read.get_tag('RZ')[0]
            if first=='A':
                self.per_cell_ta_obs[sample][ 'A_start' ] += 1
            self.per_cell_ta_obs[sample][ 'total' ] += 1

    def __repr__(self):
        return 'ScCHICLigation: no description'

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots(figsize=(4,4))

        x = []
        y = []
        for cell, cell_data in self.per_cell_ta_obs.items():
            x.append(cell_data['total'] )
            y.append( cell_data['A_start'] /  cell_data['total'] )


        ax.scatter(x,y)
        ax.set_xscale('log')
        if title is not None:
            ax.set_title(title)

        ax.set_ylabel("Fraction unique cuts starting with A")
        ax.set_xlabel("# Molecules")
        ax.set_xlim(1,None)
        ax.set_ylim(0,1)
        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()
