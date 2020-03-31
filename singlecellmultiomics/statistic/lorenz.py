import matplotlib.pyplot as plt
from singlecellmultiomics.bamProcessing import random_sample_bam
import singlecellmultiomics.pyutils as pyutils
import collections
import pandas as pd
import matplotlib
import numpy as np
import seaborn as sns
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')

def plot_lorentz(cdf, per_sample=False):
    fig, ax = plt.subplots(figsize=(6,6))
    if per_sample:
        for cell in cdf:
            ax.plot(np.linspace(0,1,cdf.shape[0]), np.cumsum(cdf[cell].fillna(0).sort_values(ascending=True))/cdf[cell].sum(),label=cell,zorder=1)

    else:
        ax.plot(np.linspace(0,1,cdf.shape[0]), np.cumsum(cdf.sum(1).fillna(0).sort_values(ascending=True))/cdf.sum().sum(),label='observed',zorder=1)
        ax.plot([0,1],[0,1],c='grey',ls=':',label='optimum',zorder=1)
        plt.title('Lorenz curve, all samples')

    ax.set_ylabel('Fraction of molecules (cumulative)')
    ax.set_xlabel('Fraction of genome')
    plt.legend()
    ax.grid(zorder=0)
    sns.despine()
    return fig, ax

class Lorenz:
    def __init__(self, args):
        pass

    def process_file(self, path):
        self.cdf = random_sample_bam(path, 10_000)

    def to_csv(self, path):
        self.cdf.to_csv(path)

    def __repr__(self):
        return f'Lorenz'

    def plot(self, target_path, title=None):
        fig, ax = plot_lorentz(self.cdf,False)

        plt.tight_layout()
        plt.savefig(target_path)
        plt.close()


        fig, ax = plt.subplots(figsize=(10,5))
        cov_per_cell = (((self.cdf>0).sum() / self.cdf.shape[0]).sort_values())
        cov_per_cell.name='fraction genome covered'
        cov_per_cell.plot.bar()

        mean_cov = cov_per_cell.mean()
        ax.axhline(mean_cov,c='green',label='mean coverage (%.3f)' % mean_cov)
        ax.set_ylabel('Fraction genome covered')
        ax.set_xlabel("Cells")
        ax.set_xticks([],[])
        sns.despine()
        plt.legend()

        plt.savefig(target_path.replace('.png', '.cell_genome_fraction.png'))
        plt.tight_layout()
        plt.close()
