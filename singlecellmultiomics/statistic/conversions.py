#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import seaborn as sns
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections
import pandas as pd
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


class ConversionMatrix(StatisticHistogram):
    def __init__(self, args, process_reads=200_000):
        StatisticHistogram.__init__(self, args)
        self.conversion_obs = collections.defaultdict(collections.Counter)
        self.base_obs = collections.defaultdict(collections.Counter)
        self.stranded_base_conversions = collections.defaultdict(
            collections.Counter)
        self.processed_reads = 0
        self.process_reads = process_reads

    def processRead(self, R1,R2):

        for read in [R1,R2]:
            if read is None:
                continue

            if self.processed_reads >= self.process_reads or read is None or read.is_unmapped or read.mapping_quality < 30:
                return

            self.processed_reads += 1
            try:
                for index, reference_pos, reference_base in read.get_aligned_pairs(
                        with_seq=True, matches_only=True):
                    query_base = read.seq[index]
                    reference_base = reference_base.upper()
                    if reference_base != 'N' and query_base != 'N':
                        k = (
                            query_base, 'R1' if read.is_read1 else (
                                'R2' if read.is_read2 else 'R?'), ('forward', 'reverse')[
                                read.is_reverse])

                        self.base_obs[reference_base][k] += 1
                        if reference_base != query_base:
                            self.conversion_obs[reference_base][k] += 1

                            if read.has_tag('RS'):
                                k = (
                                    query_base,
                                    'R1' if read.is_read1 else (
                                        'R2' if read.is_read2 else 'R?'),
                                    read.get_tag('RS'))
                                self.stranded_base_conversions[reference_base][k] += 1
            except ValueError: # Fails when the MD tag is not present
                continue

    def __repr__(self):
        return f'Observed base conversions'

    def get_df(self):
        return pd.DataFrame(
            self.base_obs).fillna(0).sort_index(
            axis=0).sort_index(
            axis=1)

    def plot(self, target_path, title=None):

        df = pd.DataFrame(
            self.conversion_obs).fillna(0).sort_index(
            axis=0).sort_index(
            axis=1)
        cm = sns.clustermap(df, row_cluster=True, col_cluster=True)
        ax = cm.ax_heatmap
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        cm.ax_heatmap.set_title('Raw conversions')
        cm.ax_heatmap.set_xlabel('reference base')
        cm.ax_heatmap.set_ylabel('sequenced base')
        plt.subplots_adjust(left=0, right=0.8)
        if title is not None:
            cm.ax_heatmap.set_title(title)
        #
        plt.savefig(target_path.replace('.png', f'.conversions.png'))
        plt.close()

        df = self.get_df()
        cm = sns.clustermap(df, row_cluster=False, col_cluster=False)
        ax = cm.ax_heatmap
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        cm.ax_heatmap.set_title('Raw conversions')
        cm.ax_heatmap.set_xlabel('reference base')
        cm.ax_heatmap.set_ylabel('sequenced base')
        #plt.tight_layout(pad=0.4, w_pad=0.8, h_pad=1.0)
        plt.subplots_adjust(left=0, right=0.8)
        if title is not None:
            cm.ax_heatmap.set_title(title)

        plt.savefig(target_path.replace('.png', f'.base_obs.png'))
        plt.close()

        if len(self.stranded_base_conversions) == 0:
            return

        df = pd.DataFrame(
            self.stranded_base_conversions).fillna(0).sort_index(
            axis=0).sort_index(
            axis=1)
        cm = sns.clustermap(df, row_cluster=False, col_cluster=False)
        ax = cm.ax_heatmap
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        cm.ax_heatmap.set_title('Raw conversions')
        cm.ax_heatmap.set_xlabel('reference base')
        cm.ax_heatmap.set_ylabel('sequenced base')
        #plt.tight_layout(pad=0.4, w_pad=0.8, h_pad=1.0)
        plt.subplots_adjust(left=0, right=0.8)
        if title is not None:
            cm.ax_heatmap.set_title(title)

        plt.savefig(
            target_path.replace(
                '.png',
                f'.RS_TAG_strand_conversions.png'))
        plt.close()

    def to_csv(self, path):

        self.get_df().to_csv(path)
