import matplotlib.patheffects as path_effects
import math
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
from collections import defaultdict, Counter
from singlecellmultiomics.utils.plotting import plot_plate
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def counter_defaultdict():
    return defaultdict(Counter)

class PlateStatistic2(object):

    def __init__(self, args):
        self.args = args

        self.stats = defaultdict(counter_defaultdict)

    def to_csv(self, path):

        export = {}
        for library, libd in self.stats.items():
            for (stat,metric), data in  libd.items():
                print(data)
                export[library,stat,metric] = data


        pd.DataFrame(
            export).to_csv(
            path.replace(
                '.csv',
                f'plate_metrics.csv'))

    def processRead(self, R1,R2=None, default_lib=None):

        for read in [R1,R2]:

            if read is None:
                continue

            cell_index = int(read.get_tag('bi'))

            if read.has_tag('MX'):
                mux= read.get_tag('MX')
                format = 384 if ('384' in mux or mux.startswith('CS2')) else 96
                if format==384:
                    n_rows = 16
                    n_cols = 24
                else:
                    n_rows = 8
                    n_cols = 12
            else:
                n_rows = 16
                n_cols = 24

            if default_lib is None:
                if read.has_tag('LY'):
                    library = read.get_tag('LY')


                else:
                    library = '_'
            else:
                library = default_lib

            cell_index = int(cell_index)-1
            row = int(cell_index/n_cols)
            col = cell_index - row*n_cols

            cell= read.get_tag('SM')
            self.stats[library][('total mapped','# reads')][(row,col, cell)] += 1

            if read.is_qcfail:
                self.stats[library]['qcfail', '# reads'][(row,col,cell)] += 1
                continue

            if read.is_duplicate:
                continue

            self.stats[library][('total mapped','# molecules')][(row,col,cell)] += 1

            if read.has_tag('lh') and  read.get_tag('lh') in ('TA','TT','AA'):
                self.stats[library][read.get_tag('lh'), 'ligated molecules'][(row,col,cell)] += 1

            break


    def __repr__(self):
        return 'Plate statistic'


    def __iter__(self):
        for library, libd in self.stats.items():
            for name, data in libd.items():
                for k,v in data.items():
                    yield (library,name,k),v


    def plot(self, target_path, title=None):
        for library, libd in self.stats.items():
            for (stat,metric), data in  libd.items():
                for log in [True,False]:

                    fig, ax, cax = plot_plate(data, log=log, cmap_name='viridis',vmin=1 if  log else 0)
                    cax.set_ylabel(metric)
                    ax.set_title(f'{stat}, {metric}')
                    plt.savefig(
                        target_path.replace(
                            '.png',
                            '') +
                        f'_{library}_{stat}_{metric.replace("# ","n_")}.{"log" if log else "linear"}.png', bbox_inches='tight', dpi=250)
                    plt.close()
            # Create oversequencing plot: (reads / molecule)
            overseq = {}
            for coordinate, n_molecules in libd[('total mapped','# molecules')].items():
                n_reads = libd[('total mapped','# reads')].get(coordinate,0)
                overseq[coordinate] = (n_reads/n_molecules) if n_reads>0 else 0

            fig, ax, cax = plot_plate(overseq, log=False, cmap_name='viridis', vmin=np.percentile(list(overseq.values()),5))
            cax.set_ylabel('reads / molecule')
            ax.set_title(f'Oversequencing')
            plt.savefig(
                target_path.replace(
                    '.png',
                    '') +
                f'_{library}_oversequencing.png', bbox_inches='tight', dpi=250)
            plt.close()

            # Create TA ratio plot: (ta molecules / total molecules)
            ta_ratio = {}
            for coordinate, n_ta_molecules in libd[('TA','ligated molecules')].items():
                n_molecules = libd[('total mapped','# molecules')].get(coordinate,0)
                ta_ratio[coordinate] = (n_ta_molecules/n_molecules) if n_molecules>0 else 0

            fig, ax, cax = plot_plate(ta_ratio, log=False, cmap_name='viridis', vmin=np.percentile(list(ta_ratio.values()),5))
            cax.set_ylabel('TA molecules / molecule')
            ax.set_title(f'TA fraction')
            plt.savefig(
                target_path.replace(
                    '.png',
                    '') +
                f'_{library}_ta_fraction.png', bbox_inches='tight', dpi=250)
            plt.close()
