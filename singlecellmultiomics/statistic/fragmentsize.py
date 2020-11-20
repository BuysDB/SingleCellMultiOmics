#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections
import seaborn as sns
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def readIsDuplicate(read):
    return (read.has_tag('RC') and read.get_tag('RC') > 1) or read.is_duplicate


class FragmentSizeHistogram(StatisticHistogram):
    def __init__(self, args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()

        self.histogramReject = collections.Counter()
        self.histogramAccept = collections.Counter()
        self.histogramMolecule = collections.Counter()

    def processRead(self, R1,R2=None):
        fragmentSize = None
        moleculeSize = None
        qcfail =False
        for read in [R1,R2]:
            if read is None:
                continue
            if read.is_qcfail:
                qcfail = True

            if read.has_tag('ms') and not read.is_duplicate:
                moleculeSize = read.get_tag('ms')

            if read.has_tag('fS'):
                fragmentSize = read.get_tag('fS')
                break

        if fragmentSize is None:
            for read in [R1,R2]:
                if read is None:
                    continue
                if read.is_qcfail:
                    qcfail = True

                if read.isize !=0:
                    fragmentSize = abs(read.isize)
                    break

        if fragmentSize is None or fragmentSize > 1_000:
            return
        #print(fragmentSize, read.reference_start,  read.reference_end,mateStart,readLen  )
        self.histogram[fragmentSize] += 1


        if qcfail:
            self.histogramReject[fragmentSize] += 1

        else:
            if moleculeSize is not None:
                self.histogramMolecule[fragmentSize] += 1
            self.histogramAccept[fragmentSize] += 1


    def __repr__(self):
        return f'The average fragment size is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'

    def plot(self, target_path, title=None):
        fig, ax = plt.subplots()
        ax.bar(
            list(
                self.histogram.keys()), list(
                self.histogram.values()), width=1)
        if title is not None:
            ax.set_title(title)

        ax.set_xlabel("Fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        plt.savefig(target_path)
        sns.despine()
        plt.close()

        fig, ax = plt.subplots()
        plt.bar(
            list(
                self.histogramReject.keys()), list(
                self.histogramReject.values()), width=1, color='r')
        if title is not None:
            plt.title(title)

        ax.set_xlabel("Rejected fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        sns.despine()
        plt.savefig(target_path.replace('.png', '.rejected.png'))
        plt.close()

        fig, ax = plt.subplots()
        ax.bar(
            list(
                self.histogramAccept.keys()), list(
                self.histogramAccept.values()), width=1, color='g')
        if title is not None:
            plt.title(title)


        ax.set_xlabel("Accepted fragment size [bp]")
        ax.set_ylabel("# Fragments")
        plt.tight_layout()
        sns.despine()
        plt.savefig(target_path.replace('.png', '.accepted.png'))
        plt.close()


        fig, ax = plt.subplots()
        ax.bar(
            list(
                self.histogramMolecule.keys()), list(
                self.histogramMolecule.values()), width=1, color='g')
        if title is not None:
            plt.title(title)

        sns.despine()
        ax.set_xlabel("Molecule size [bp]")
        ax.set_ylabel("# Molecules")
        plt.tight_layout()
        plt.savefig(target_path.replace('.png', '.molecules.png'))
        plt.close()
