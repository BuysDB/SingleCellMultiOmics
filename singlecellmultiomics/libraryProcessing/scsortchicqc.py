#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from singlecellmultiomics.utils import createRowColorDataFrame
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import squareform
import scipy.stats
import pandas as pd
import numpy as np
import glob
import seaborn as sns
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import json
from singlecellmultiomics.utils.plotting import plot_plate, plot_plate_layout
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
import os
from singlecellmultiomics.libraryProcessing import SampleSheet

def read_plate_statistics(path):

    stat_tab = defaultdict(dict)
    for (_lib,statistic,(x,y,cell)),val in pd.Series( pd.read_pickle(path)['PlateStatistic2'] ).items():
        stat_tab[statistic][cell] = val

    stat_tab = pd.DataFrame( stat_tab ).fillna(0)
    return stat_tab

def read_count_table(path):
    library = path.split('/')[-3]
    df = pd.read_pickle(path).T
    # rebuild index
    df.index = [ f'{library}_{int(cell)+1}' for cell in df.index]
    return df

def sample_sheet_to_condition_labels(sample_sheet):
    cell_labels = []
    layout_formats = sample_sheet['layout_format']
    for library, layout_name in sample_sheet["library_layout"].items():
        #print(well2coord[layout_formats[layout_name]])
        #well2index= {k:tuple(v) for k,v in well2coord[layout_formats[layout_name]].items()}

        layout_format = layout_formats[layout_name]

        # Annotate the cells by well information:
        cell_labels.append(
            pd.Series(
            {
                f'{library}_{sample_sheet["well2index"][layout_format][well]}': well_label
                for well, well_label in sample_sheet["layouts"][layout_name].items()
                }
            )
        )
    return pd.concat(cell_labels)

def read_contaminant_info(sortchicstats_paths):
    # Read contaminant information:
    contaminant_info = []
    for p in sortchicstats_paths:
        with open(p) as f:
            ci = json.load(f)
            contaminant_info.append(
                (pd.DataFrame(ci['sc_se_scaffold_cov']).T / pd.Series( ci['sc_se_reads'])).T.fillna(0)
            )
    contaminant_info = pd.concat(contaminant_info)
    contaminant_info.columns = pd.MultiIndex.from_tuples( [(col,'fraction of reads') for col in contaminant_info] )
    return contaminant_info

if __name__ == '__main__':

    import matplotlib
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['font.family'] = 'Helvetica'

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract cut distribution from bam file')

    argparser.add_argument('sample_sheet', type=str)
    argparser.add_argument('count_tables_sortchicstats_statistics', type=str, nargs='+')
    argparser.add_argument('-o', default='quality_control',type=str)
    argparser.add_argument('-ignore_marks',type=str, help='Comma separated marks to ignore')
    #argparser.add_argument('--per_mark', action='store_true',help='Perform quality control predictions per mark, requires enough cells being available, one plate is often not sufficient')
    args = argparser.parse_args()

    ignore_marks = None if args.ignore_marks is None else args.ignore_marks.split(',')

    # decompose count_tables_sortchicstats_statistics
    sortchicstats_paths = []
    statistics_paths = []
    count_table_paths = []
    for path in args.count_tables_sortchicstats_statistics:

        if path.endswith('statistics.pickle.gz'):
            statistics_paths.append(path)
        elif path.endswith('sortchicstats.json'):
            sortchicstats_paths.append(path)
        elif path.endswith('.pickle.gz'):
            count_table_paths.append(path)
        else:
            raise ValueError(f'File format of {path} not understood')


    sample_sheet = SampleSheet(args.sample_sheet)

    # Read the count tables
    df = pd.concat([read_count_table(path) for path in count_table_paths])
    # Add mark as first level of df, library second, cell third
    df.index = pd.MultiIndex.from_tuples([(sample_sheet['marks'][cell.split('_')[0]], cell.split('_')[0], int(cell.split('_')[1])) for cell in df.index])

    avail_marks = df.index.get_level_values(0).unique()
    print('Target marks:')
    print(avail_marks)

    if ignore_marks is not None:
        print('Ignoring', ignore_marks)
        df = df.drop(ignore_marks,level=0)
        sample_sheet.drop_mark(ignore_marks)
        avail_marks = df.index.get_level_values(0).unique()
        print('Remaining target marks:')
        print(avail_marks)


    # Find for each plate, for each cell the cell idx -> the class
    cell_labels = sample_sheet_to_condition_labels(sample_sheet)
    contaminant_info = read_contaminant_info(sortchicstats_paths)

    # Perform classification:
    plate_stats = pd.concat([read_plate_statistics(p) for p in statistics_paths])
    y = cell_labels=='empty'
    rf = RandomForestClassifier(class_weight='balanced')

    X = plate_stats.loc[y.index]
    X[('AA', 'ligated molecules')]/=X[('total mapped',       '# molecules')]
    X[('TA', 'fraction ligated molecules')]= X[('TA', 'ligated molecules')] / X[('total mapped',       '# molecules')]
    X[('TT', 'ligated molecules')]/=X[('total mapped',       '# molecules')]
    X[('qcfail', '# reads')]/=X[('total mapped',       '# molecules')]
    X[('duprate', 'pct')] =X[('total mapped',       '# molecules')]/X[('total mapped', '# reads')]

    y[X[('total mapped','# reads')]<500] = True
    X = X.join(contaminant_info)

    predictions = []
    for train_index, test_index in KFold(n_splits=8, shuffle=True, random_state=None).split(X):
        rf.fit(X.iloc[train_index],y.iloc[train_index])
        predictions.append( pd.Series(rf.predict_proba(X.iloc[test_index])[:,0],index=X.iloc[test_index].index))
    predictions = pd.concat(predictions)

    layout_format = 'scChIC384'
    index2well = {v:k for k,v in sample_sheet["well2index"][layout_format].items()}

    library_coord_verdict = defaultdict(dict)
    library_coord_posterior = defaultdict(dict)

    # find threshold such that all empty wells are classified as empty.
    empty_posteriors = []
    for cell, pred in predictions.items():
        if y[cell]:
            empty_posteriors.append(pred)
    th = np.percentile(empty_posteriors,99)

    qc_pass =  {}

    for cell, pred in predictions.items():
        library, cidx = cell.split('_')
        cidx=int(cidx)

        cell_qc_passed = pred>=th if not y[cell] else False
        qc_pass[cell] = cell_qc_passed
        library_coord_verdict[library][tuple(sample_sheet['well2coord']['scChIC384'][index2well[cidx]])] = cell_qc_passed
        library_coord_posterior[library][tuple(sample_sheet['well2coord']['scChIC384'][index2well[cidx]])] = pred if not y[cell] else 0
    qc_pass = pd.Series(qc_pass)

    ## Create tables:
    table_folder = f'{args.o}/tables'
    if not os.path.exists(table_folder):
        os.makedirs(table_folder)

    qc_pass.name = 'cell_is_qc_pass'
    qc_pass.to_csv(f'{table_folder}/qc_pass.csv')

    ## Create plots:
    plotpath = f'{args.o}/plots/QC_scatter'
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    for xvar,yvar in [
    [('TA', 'fraction ligated molecules',"linear"), ('total mapped', '# molecules',"log")],
    [('total mapped', '# reads',"log"), ('total mapped', '# molecules',"log")],
     [('E. coli RHB09-C15','fraction of reads',"log"), ('total mapped', '# molecules',"log")]

    ]:

        fig,ax = plt.subplots(figsize=(4,3))
        mpb = ax.scatter(X[~y][xvar[:2]],X[~y][yvar[:2]],c=predictions.loc[X[~y].index],alpha=0.9,label='well')
        # Show locations of empty wells:
        ax.scatter(X[y][xvar[:2]],X[y][yvar[:2]],c='k',marker='x',s=50,label='empty well')
        ax.set_yscale(yvar[2])
        ax.set_xscale(xvar[2])
        plt.legend()
        plt.xlabel(f'{xvar[0]} [{xvar[1]}]')
        plt.ylabel(f'{yvar[0]} [{yvar[1]}]')
        #plt.ylabel('total mapped molecules')
        axcol = plt.colorbar(mpb)
        axcol.set_label('Quality score')


        descriptor = ('-'.join(xvar) + '_vs_' + '-'.join(yvar)).replace('#','').replace(' ','_')
        plt.tight_layout()
        plt.savefig( f'{plotpath}/{descriptor}.QC_scatter.png' )
        plt.savefig( f'{plotpath}/{descriptor}.QC_scatter.svg' )


    plotpath = f'{args.o}/plots/QC_plate_score'
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    for lib, d in library_coord_posterior.items():
        fig,ax,cbar = plot_plate(d,vmax=1,vmin=0,usenorm=True,log=False, cmap_name='viridis') #,suptitle=f'{lib}')
        cbar.set_position([0.98,0.2,0.05,0.1])
        cbar.set_title('Quality score')
        st = plt.suptitle(lib)
        #plt.tight_layout()

        fig.subplots_adjust(wspace=0.05)
        plt.savefig( f'{plotpath}/{lib}.QC_plate_score.png',bbox_extra_artists=[st,cbar],bbox_inches='tight')
        plt.savefig( f'{plotpath}/{lib}.QC_plate_score.svg',bbox_extra_artists=[st,cbar],bbox_inches='tight')



    plotpath = f'{args.o}/plots/QC_plate'
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    for lib, d in library_coord_verdict.items():
        fig,ax,cbar = plot_plate(d,vmax=1,vmin=0,usenorm=True,log=False, cmap_name='viridis') #,suptitle=f'{lib}')
        cbar.set_position([0.98,0.2,0.05,0.1])
        cbar.set_title('QC pass')
        st = plt.suptitle(lib)
        #plt.tight_layout()

        fig.subplots_adjust(wspace=0.05)
        plt.savefig( f'{plotpath}/{lib}.QC_plate.png',bbox_extra_artists=[st,cbar],bbox_inches='tight')
        plt.savefig( f'{plotpath}/{lib}.QC_plate.svg',bbox_extra_artists=[st,cbar],bbox_inches='tight')

    # Perform correlation analysis
