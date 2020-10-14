#!/usr/bin/env python
# -*- coding: utf-8 -*-
from multiprocessing import Pool
import pysam
import pandas as pd
import argparse
from collections import Counter, defaultdict
import numpy as np
import seaborn as sns
import os
import gzip
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
mpl.rcParams['figure.dpi'] = 300

if __name__=='__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Plot base calling covariates""")
    argparser.add_argument('covariates_in', metavar='covariates_in', type=str, nargs='+')
    argparser.add_argument('-o', help='output folder',  type=str, default='./plots')

    args = argparser.parse_args()



    if not os.path.exists(args.o):
        os.makedirs(args.o)

    for path in args.covariates_in:
        alias = os.path.basename(path).replace('.pickle.gz','')
        print(f'Loading covariates from {path} ')
        with gzip.open(path,'rb') as i:
            covariates = pickle.load(i)

        # Integrate over cycle
        for mate_select in (True,False):
            try:
                qf = pd.DataFrame([(cycle, neg/(pos+neg), context) for (qual, mate, cycle, context), (pos, neg) in covariates.items() if mate==mate_select and qual>=30])
                qf.columns=['Cycle','Accuracy','Context']
                qf = qf.sort_values('Cycle')
                sns.boxplot(data=qf,x=qf.columns[0],y=qf.columns[1])
                plt.savefig(f'{args.o}/cycle_r{"R2" if mate_select else "R1"}_{alias}.png')
                plt.close()
            except Exception as e:
                pass
        # Integrate over contexts:
        context_acc = pd.DataFrame([ (c,v['Accuracy'].mean())  for c,v in qf.groupby('Context') if len(c)==3]).sort_values(1).set_index(0)
        context_acc = context_acc.sort_values(1)

        fig, axes = plt.subplots(2,1, figsize=(16,10))

        ax = axes[0]
        for c in context_acc.index[:10]:
            v = qf[qf['Context']==c]

            v.groupby('Cycle').mean().plot.line(y='Accuracy',label=c, ax=ax) #x='Cycle',
        ax.set_title('10 worst performing 3bp contexts\nfor base calls with high phred score')
        sns.despine()
        ax.grid()

        ax = axes[1]
        for c in context_acc.index[-10:]:
            v = qf[qf['Context']==c]
            v.groupby('Cycle').mean().plot.line(y='Accuracy',label=c, ax=ax) #x='Cycle',
        ax.set_title('10 best performing 3bp contexts\nfor base calls with high phred score')
        sns.despine()
        ax.grid()
        plt.savefig(f'{args.o}/contexts_cycle_{alias}.png')
        plt.close()

        # Zoomins:
        for mate_select in (True, False):
            try:
                qf = pd.DataFrame(
                    [(cycle, neg / (pos + neg), context) for (qual, mate, cycle, context), (pos, neg) in covariates.items()
                     if mate == mate_select and qual >= 36 and cycle < 65])
                qf.columns = ['Cycle', 'Accuracy', 'Context']
                qf = qf.sort_values('Cycle')
                sns.boxplot(data=qf, x=qf.columns[0], y=qf.columns[1])
                sns.despine()
                ax.grid()

                plt.savefig(f'{args.o}/contexts_cycle_zoomed_{alias}.png')
                plt.ylim(0.996, 1.0001)

                plt.show()
            except Exception as e:
                pass
