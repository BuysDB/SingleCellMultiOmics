#!/usr/bin/env python
# -*- coding: utf-8 -*-
from multiprocessing import Pool
import pysam
import pandas as pd
import os
from scipy.optimize import curve_fit
import argparse
from singlecellmultiomics.bamProcessing.bamFunctions import get_contigs_with_reads, get_r1_counts_per_cell
from collections import Counter, defaultdict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class DivCounter(Counter):
    """Divisable counter"""
    def __truediv__(self,other):
        result = Counter()
        for k,v in self.items():
            result[k] = v/other
        return result


def _get_sc_cut_dictionary(args):

    bam, contig, strand_specific, filter_function = args
    cut_positions = defaultdict(Counter)
    with pysam.AlignmentFile(bam) as alignments:
        for read in alignments.fetch(contig):

            if not filter_function(read):
                continue
            cut_positions[read.get_tag('SM')][
                (read.is_reverse, read.get_tag('DS'))
                    if strand_specific else
                read.get_tag('DS')
                ]+=1

    return contig,cut_positions


def read_counts_function(read):
    if not read.is_read1 or read.is_duplicate or read.is_qcfail or read.mapping_quality==0:
        return False
    return True

def get_sc_cut_dictionary(bam_path: str, filter_function=None, strand_specific=False):
    """
    Generates cut distribution dictionary  (contig)->sample->position->obs

    """
    if filter_function is None:
        filter_function = read_counts_function
    cut_sites = {}
    with Pool() as workers:
        with pysam.AlignmentFile(bam_path) as alignments:
            for contig,r in workers.imap_unordered(
                _get_sc_cut_dictionary, (
                    (bam_path, contig,strand_specific, filter_function)
                    for contig in get_contigs_with_reads(bam_path) )):
                cut_sites[contig]=r

    return cut_sites


def dictionary_to_diff_vector(d,sample: str, vmin: float, vmax: float):
    """Convert a dict {contig:sample:position:obs} into sorted vector [ distance, distance, ..]"""
    return np.array([
        v for v in np.clip(
        np.concatenate(
            [np.diff(sorted(d[contig][sample])) for contig in d])
            ,vmin,vmax) if v>vmin and v<vmax])


def analyse(bam_path,output_dir, create_plot=False, min_distance=20, max_distance=800, verbose=False, strand_specific=False):

    if verbose:
        print('Obtaining molecules per cell .. ', end='\r')
    cpr = get_r1_counts_per_cell(bam_path)

    if verbose:
        print('Molecules per cell:            ')
        for cell, obs in cpr.most_common():
            print(f'\t{cell}\t{obs}')

    if verbose:
        print('Obtaining cuts per cell .. ', end='\r')

    cut_sites = get_sc_cut_dictionary(bam_path, strand_specific=strand_specific)


    all_counts = {}
    for cell, total_molecules in cpr.most_common():
        # Write from 0 to max_distance table
        all_counts[cell] = DivCounter(dictionary_to_diff_vector(cut_sites,cell,0,max_distance))

    cut_count_df = pd.DataFrame(all_counts).sort_index().sort_index(1).fillna(0)
    cut_count_df.to_csv(f'{output_dir}/counts.csv')

    if verbose:
        print('Obtaining cuts per cell [ OK ]')
        print('Fitting and plotting ..', end='\r')

    if create_plot:
        try:
            cut_count_df.index.name='distance between cuts'
            filtered_count_df = cut_count_df.loc[:, cut_count_df.sum()>100]
            sns.clustermap((filtered_count_df / filtered_count_df.loc[20:].mean()).T,
                           cmap='viridis', vmax=3,
                           metric='correlation', col_cluster=False,
                           method='ward',figsize=(8,20))
            plt.tight_layout()
            plt.savefig(f'{output_dir}/heatmap.png')


            #ax.figure.subplots_adjust(left=0.3)  # change 0.3 to suit your needs.

        except Exception as e:
            print(e)


    def function_to_fit(xdata, period, offset, amplitude, decay, mean ):
        frequency = 1/period
        return  (amplitude*np.cos((2*np.pi*(frequency)*(xdata+offset) ))) * np.exp(-xdata*(1/decay)) + mean

    # Bounds for fitting:
    bounds=(
        (150,300), # Frequency (b)
        (-30,30), # offset (c)
        (1,400), # amplitude
        (100,1900), # decay
        (1,99999), # mean
    )

    if create_plot:
        sc_plot_dir = f'{output_dir}/sc_plots'
        if not os.path.exists(sc_plot_dir):
            os.makedirs(sc_plot_dir)

    smooth_small_signals = {}
    smooth_big_signals = {}
    fit_params_per_cell = defaultdict(dict)
    for cell, total_molecules in cpr.most_common():
        try:
            sc_counts = pd.DataFrame({
                cell:DivCounter(
                    dictionary_to_diff_vector(cut_sites,cell,min_distance,max_distance))})

            if create_plot:
                fig, ax = plt.subplots(figsize=(10,3))

            big_window = 35
            smooth = sc_counts.rolling(window=big_window,center=True).mean()
            smooth_big_signals[cell] = smooth[cell]
            if create_plot:
                ax.plot(smooth.index, smooth[cell],label=f'{big_window}bp sliding window')
                limits = ax.get_ylim()

            xdata = sc_counts[cell].index
            ydata = sc_counts[cell].values
            if len(ydata)==0:
                continue

            xdata = xdata[~np.isnan(ydata)]
            ydata = ydata[~np.isnan(ydata)]
            fit_params = curve_fit(function_to_fit, xdata, ydata,bounds=(np.array(bounds).T[0], np.array(bounds).T[1]))[0]

            if create_plot:
                plt.scatter(xdata,ydata, c='grey', s=1, label='Raw data')

            period, offset, amplitude, decay,mean = fit_params
            fit_params_per_cell['period'][cell] = period
            fit_params_per_cell['offset'][cell] = offset
            fit_params_per_cell['amplitude'][cell]= amplitude
            fit_params_per_cell['decay'][cell] = decay
            fit_params_per_cell['mean'][cell] = mean

            if not create_plot:
                continue

            plt.plot(xdata,function_to_fit(xdata,*fit_params), c='r',
                     label=f'Fit : per:{period:.0f} ph:{offset:.0f} mean:{mean:.0f} dec:{decay:.2f}')

            ax.axhline(mean,c='k')
            ax.axvline(period-offset,c='b',lw=1)
            ax.axvline(2*period-offset,c='b',lw=1)

            ax.set_title(f'{cell},\n{total_molecules} molecules' )
            ax.set_xlabel(f'distance to nearest cut (bp)' )
            ax.set_ylabel(f'# cuts' )

            ax.set_ylim( (limits[0]*0.9,limits[1]*1.1))

            sns.despine()
            ax.grid()
            plt.legend()

            plt.tight_layout()

            plt.savefig(f'{sc_plot_dir}/{cell}.png')
            plt.close()

            # Plot residual with smoothed function


        except RuntimeError as e:
            print(f'Could not fit data for {cell}, ( {total_molecules} molecules )')
            pass

    if verbose:
        print('Fitting and plotting [ OK ]')
        print('Writing files ..', end='\r')

    # Write tables to disk
    tmp = {'molecules_total':cpr}
    tmp.update(fit_params_per_cell)
    df = pd.DataFrame(tmp)
    df.to_csv(f'{output_dir}/fit.csv')


    if verbose:
        print('All done     ')

if __name__ == '__main__':

    import matplotlib
    matplotlib.rcParams['figure.dpi'] = 160
    matplotlib.use('Agg')

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract cut distribution from bam file')

    argparser.add_argument('alignmentfile', type=str)
    argparser.add_argument('-o', type=str, required=True, help='Output folder')
    argparser.add_argument('-max_distance', type=int,default=800, help='Maximum distance in both plots and output tables')
    args = argparser.parse_args()

    if not os.path.exists(args.o):
        os.makedirs(args.o)

    analyse(args.alignmentfile, args.o, create_plot=True, verbose=True,strand_specific=False,max_distance=args.max_distance)
