#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections
import pandas as pd

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


class RejectionReasonHistogram(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

    def processRead(self, R1,R2):

        for read in [R1,R2]:
            if read is None:
                continue
            if read.has_tag('RR'):
                self.histogram[read.get_tag('RR')] += 1
            break

    def __repr__(self):
        rt = 'Rejection reasons:'
        for reason, obs in self.histogram.most_common():
            rt += f'{reason}\t:\t{obs}\n'
        return rt

    def __iter__(self):
        return iter(self.histogram.most_common())

    def plot(self, target_path, title=None):
        if len(self.histogram) == 0:
            print('Not plotting rejection reasons, rejection tag RR was never seen')
            return

        df = pd.DataFrame.from_dict({'Reason': dict(self)}).T

        df.plot.bar(figsize=(10, 4)).legend(bbox_to_anchor=(1, 0.98))
        if title is not None:
            plt.title(title)

        plt.tight_layout()
        plt.subplots_adjust(right=0.6)
        plt.savefig(target_path)

        ax = plt.gca()
        ax.set_yscale('log')
        plt.savefig(target_path.replace('.png', '.log.png'))
        plt.close()
