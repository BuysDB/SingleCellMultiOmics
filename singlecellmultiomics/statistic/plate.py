import matplotlib.patheffects as path_effects
import math
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def human_readable(value, targetDigits=2, fp=0):

    # Float:
    if value < 1 and value > 0:
        return('%.2f' % value)

    if value == 0.0:
        return('0')

    baseId = int(math.floor(math.log10(float(value)) / 3.0))
    suffix = ""
    if baseId == 0:
        sVal = str(round(value, targetDigits))
        if len(sVal) > targetDigits and sVal.find('.'):
            sVal = sVal.split('.')[0]

    elif baseId > 0:

        sStrD = max(0, targetDigits -
                    len(str('{:.0f}'.format((value / (math.pow(10, baseId * 3)))))))

        sVal = ('{:.%sf}' % min(fp, sStrD)).format(
            (value / (math.pow(10, baseId * 3))))
        suffix = 'kMGTYZ'[baseId - 1]
    else:

        sStrD = max(0, targetDigits -
                    len(str('{:.0f}'.format((value * (math.pow(10, -baseId * 3)))))))
        sVal = ('{:.%sf}' % min(fp, sStrD)).format(
            (value * (math.pow(10, -baseId * 3))))
        suffix = 'mnpf'[-baseId - 1]

        if len(sVal) + 1 > targetDigits:
            # :(
            sVal = str(round(value, fp))[1:]
            suffix = ''

    return('%s%s' % (sVal, suffix))


# Visualize the following:
# PER LIBRARY / DEMUX method
# total fragments
# total fragments with correct site
# unique molecules

# 384 well format:

well2index = collections.defaultdict(dict)
index2well = collections.defaultdict(dict)
rows = string.ascii_uppercase[:16]
columns = list(range(1, 25))

for ci in range(1, 385):
    i = ci - 1
    rowIndex = math.floor(i / len(columns))
    row = rows[rowIndex]
    column = columns[i % len(columns)]
    well2index[384][(row, column)] = ci
    index2well[384][ci] = (row, column)


rows96 = string.ascii_uppercase[:8]

columns96 = list(range(1, 13))

for ci in range(1, 97):
    i = ci - 1
    rowIndex = math.floor(i / len(columns96))
    row = rows96[rowIndex]
    column = columns96[i % len(columns96)]
    well2index[96][(row, column)] = ci
    index2well[96][ci] = (row, column)


class PlateStatistic(object):

    def __init__(self, args):
        self.args = args

        self.rawFragmentCount = collections.defaultdict(
            collections.Counter)  # (library, mux) -> cell -> counts
        self.usableCount = collections.defaultdict(
            collections.Counter)  # (library, mux) -> cell -> counts
        self.moleculeCount = collections.defaultdict(
            collections.Counter)  # (library, mux) -> cell -> counts
        self.skipReasons = collections.Counter()

    def to_csv(self, path):
        pd.DataFrame(
            self.moleculeCount).to_csv(
            path.replace(
                '.csv',
                'molecules.csv'))
        pd.DataFrame(
            self.usableCount).to_csv(
            path.replace(
                '.csv',
                'usable_reads.csv'))
        pd.DataFrame(
            self.rawFragmentCount).to_csv(
            path.replace(
                '.csv',
                'raw_fragments.csv'))

    def processRead(self, R1,R2):

        for read in [R1,R2]:

            if read is None:
                continue

            if not read.has_tag('MX'):
                return

            self.rawFragmentCount[(read.get_tag('LY'),
                                   read.get_tag('MX'))][read.get_tag('SM')] += 1

            if read.get_tag('MX').startswith('CS2'):
                if read.has_tag('XT') or read.has_tag('EX'):
                    if read.is_read1:  # We only count reads2
                        return
                    self.usableCount[(read.get_tag('LY'),
                                      read.get_tag('MX'))][read.get_tag('SM')] += 1

                    if read.has_tag('RC') and read.get_tag('RC') == 1:
                        self.moleculeCount[(read.get_tag('LY'), read.get_tag(
                            'MX'))][read.get_tag('SM')] += 1
            else:

                if read.has_tag('DS'):
                    if not read.is_read1:
                        self.skipReasons['Not R1'] += 1
                        return

                    self.usableCount[(read.get_tag('LY'),
                                      read.get_tag('MX'))][read.get_tag('SM')] += 1
                    if not read.is_duplicate:
                        self.moleculeCount[(read.get_tag('LY'), read.get_tag(
                            'MX'))][read.get_tag('SM')] += 1
                else:
                    self.skipReasons['No DS'] += 1
            break


    def __repr__(self):
        return 'Plate statistic'

    def cell_counts_to_dataframe(self, cell_counts, mux, name='raw_reads'):
        df = pd.DataFrame({name: cell_counts})

        offset = 0 # Offset is zero for all protocols since 0.1.12

        format = 384 if ('384' in mux or mux.startswith('CS2')) else 96

        df['col'] = [index2well[format]
                     [(offset + int(x.rsplit('_')[-1]))][1] for x in df.index]
        df['row'] = [-rows.index(index2well[format]
                                 [(offset + int(x.rsplit('_')[-1]))][0]) for x in df.index]
        df['size'] = (df[name] / np.percentile(df[name], 99) * 200)

        return df

    def __iter__(self):
        for data, name in [
            (self.rawFragmentCount, 'raw_reads'),
            (self.usableCount, 'usable_reads'),
                (self.moleculeCount, 'unique_molecules')]:
            for (library, mux), cellCounts in data.items():
                df = self.cell_counts_to_dataframe(cellCounts, mux, name=name)
                for i, row in df.iterrows():
                    yield i, row

    def plot(self, target_path, title=None):
        for data, name in [
            (self.rawFragmentCount, 'raw_reads'),
            (self.usableCount, 'usable_reads'),
                (self.moleculeCount, 'unique_molecules')]:

            for (library, mux), cellCounts in data.items():

                df = self.cell_counts_to_dataframe(cellCounts, mux, name=name)
                df.plot.scatter(x='col', y='row', s=df['size'],
                                c=[(0.2, 0.2, 0.5, 0.9)]
                                )

                # Annotate the outliers with values:
                ax = plt.gca()
                for ii, row in df.iterrows():
                    if row[name] > 0 and (
                        row[name] < np.percentile(
                            df[name],
                            5) or row[name] > np.percentile(
                            df[name],
                            95)):
                        text = ax.annotate(human_readable(int(row[name])), (row['col'], row['row']),
                                           ha='center', va='baseline', color='w', size=7)
                        text.set_path_effects([path_effects.Stroke(
                            linewidth=3, foreground='black'), path_effects.Normal()])

                plt.yticks(
                    sorted(
                        df['row'].unique())[
                        ::-1],
                    sorted(rows),
                    rotation=0)
                plt.xticks(
                    sorted(
                        df['col'].unique()),
                    sorted(columns),
                    rotation=0)
                plt.title(fr'{name} with ${mux}$ adapter' + f'\n{library}')

                # Create legend:
                #ld = []
                # for x in np.linspace(1, max(df[name]), 4):
            #        size = (x/np.percentile(df[name],99))*200
                #    ld.append( mlines.Line2D([], [], color='blue', marker='.', linestyle='None',
                # markersize=np.sqrt(size), label=f'{int(x)}:{size}'))
                # plt.legend(handles=ld)

                plt.savefig(
                    target_path.replace(
                        '.png',
                        '') +
                    f'{name}_{mux}_{library}.png')
                plt.close()
