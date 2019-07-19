#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections
import pandas as pd
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class MethylationContextHistogram(StatisticHistogram):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histograms_three = collections.Counter()
        self.histograms_pent = collections.Counter()
        self.histograms_three_un = collections.Counter()
        self.histograms_pent_un = collections.Counter()


    def processRead(self,read):
        if not read.is_paired or not read.is_read1 or read.is_duplicate:
            return

        if read.has_tag('Qm'):
            for modified_context in read.get_tag('Qm').split(','):
                self.histograms_pent[modified_context]+=1

        if read.has_tag('Qu'):
            for modified_context in read.get_tag('Qu').split(','):
                self.histograms_pent_un[modified_context]+=1

        if read.has_tag('Cu'):
            for modified_context in read.get_tag('Cu').split(','):
                self.histograms_three_un[modified_context]+=1

        if read.has_tag('Cm'):
            for modified_context in read.get_tag('Cm').split(','):
                self.histograms_three[modified_context]+=1

    def __repr__(self):
        return f'Methylation status.'

    def plot(self, target_path, title=None):



        for d, name, ncol in [(self.histograms_three_un,'3bp_context_unmodified',1),
            (self.histograms_three,'3bp_context_modified',1),
            (self.histograms_pent_un,'5bp_context_unmodified',3),
            (self.histograms_pent,'5bp_context_modified',3)]:
            if len(d)<1:
                print(f'No methylation data [{name}], not making plot')
                continue
            df = pd.DataFrame.from_dict({'obs':d}).T
            df.plot.bar(figsize=(15,6)).legend(bbox_to_anchor=(1, 0.98), ncol=3)
            if title is not None:
                plt.title(title)

            plt.tight_layout()
            plt.subplots_adjust(right=0.6)
            plt.savefig(target_path.replace('.png',f'.{name}.png'))

            ax = plt.gca()
            ax.set_yscale('log')
            plt.savefig(target_path.replace('.png',f'{name}.log.png'))
            plt.close()


    def to_csv(self, path):

        pd.DataFrame({'obs':self.histograms_three_un}).to_csv(path.replace('.csv','unmodified_3_base_context.csv'))
        pd.DataFrame({'obs':self.histograms_pent_un}).to_csv(path.replace('.csv','unmodified_5_base_context.csv'))
        pd.DataFrame({'obs':self.histograms_three}).to_csv(path.replace('.csv','modified_3_base_context.csv'))
        pd.DataFrame({'obs':self.histograms_pent}).to_csv(path.replace('.csv','modified_5_base_context.csv'))
