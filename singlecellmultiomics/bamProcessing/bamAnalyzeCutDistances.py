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
import seaborn as sns
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

class DivCounter(Counter):
    """Divisable counter"""
    def __truediv__(self,other):
        result = Counter()
        for k,v in self.items():
            result[k] = v/other
        return result

def find_nearest(array, values):
    idxes = np.searchsorted(array, values, side="left")
    r = []
    for value, idx in zip(values, idxes):
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
            r.append(array[idx - 1])
        else:
            r.append(array[idx])
    return r

def calculate_distance(vector_target: np.array, vector_viewpoint: np.array, max_range: float):
    # Calculate distance between viewpoint and target, skip locations with a nan (will not be returned in the result)
    existing = ~(np.isnan(vector_viewpoint) & ~np.isnan(vector_target))
    if existing.sum() == 0:
        return []
    dist = vector_viewpoint[existing] - vector_target[existing]
    return dist[(dist > -max_range) * (dist < max_range)]

def dictionary_to_diff_vector(d,sample: str, vmin: float, vmax: float):
    """Convert a dict {contig:sample:position:obs} into sorted vector [ distance, distance, ..]"""
    return np.array([
        v for v in np.clip(
        np.concatenate(
            [np.diff(sorted(d[contig][sample])) for contig in d])
            ,vmin,vmax) if v>vmin and v<vmax])

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

def extract_indices(haystack, indices, fill):
    return np.array([haystack[index] if index > 0 and index < len(haystack) else np.nan for index in indices])


def find_nearest_above(needles, haystack):
    indices = np.searchsorted(haystack, needles, side="right")
    return extract_indices(haystack, indices, np.nan)


def find_nearest_below(needles, haystack):
    haystack_rev = -haystack
    haystack_rev.sort()
    indices = np.searchsorted(haystack_rev, -needles, side="right")
    return np.abs(extract_indices(haystack_rev, indices, np.nan))


def get_stranded_pairwise_counts(sc_cut_dict_stranded, max_range=3000):
    """
    Obtain how many observations exist of different types of pairs of molecules

    Args:
        sc_cut_dict_stranded(dict) : { contig: { sample: { Counter( position: obs ) .. }}}
        max_range(int) : maximum distance to record

    Returns:
            distance_counter_fwd_above
         distance_counter_fwd_below
            distance_counter_rev_above
            distance_counter_rev_below
    """
    distance_counter_fwd_above = defaultdict(Counter)
    distance_counter_fwd_below = defaultdict(Counter)
    distance_counter_rev_above = defaultdict(Counter)
    distance_counter_rev_below = defaultdict(Counter)

    for contig in sc_cut_dict_stranded:

        for sample in sc_cut_dict_stranded[contig].keys():
            forward = np.array([pos for strand, pos in sc_cut_dict_stranded[contig][sample] if not strand])
            reverse = np.array([pos for strand, pos in sc_cut_dict_stranded[contig][sample] if strand])

            if len(forward) <= 1 or len(reverse) <= 1:
                continue

            forward.sort()
            reverse.sort()

            # for each position on the fwd strand find the closest fragment on the forward strand.
            # [>>>>>>>> .....|
            #         <<<<<<<

            nearest_fwd_above = find_nearest_above(forward, reverse)
            distance_counter_fwd_above[sample] += Counter(calculate_distance(forward, nearest_fwd_above, max_range))

            #              >>>>>>>>
            #   <<<<<<<

            nearest_fwd_below = find_nearest_below(forward, reverse)
            distance_counter_fwd_below[sample] += Counter(calculate_distance(forward, nearest_fwd_below, max_range))

            #  >>>>>>> ..........|
            #              <<<<<<]
            nearest_rev_above = find_nearest_above(reverse, forward)
            distance_counter_rev_above[sample] += Counter(calculate_distance(reverse, nearest_rev_above, max_range))

            #              >>>>>>>>
            #   <<<<<<<

            nearest_rev_below = find_nearest_below(reverse, forward)
            distance_counter_rev_below[sample] += Counter(calculate_distance(reverse, nearest_rev_below, max_range))

    return distance_counter_fwd_above, distance_counter_fwd_below, distance_counter_rev_above, distance_counter_rev_below

def read_counts_function(read):
    if not read.is_read1 or read.is_duplicate or read.is_qcfail or read.mapping_quality==0:
        return False
    return True


def strict_read_counts_function(read):
    if not read.is_read1 or \
        read.is_duplicate or \
        read.is_qcfail or \
        read.mapping_quality<50 or \
        'S' in read.cigarstring or \
        'I' in read.cigarstring or \
        not read.is_proper_pair or \
        read.get_tag('NM')>1:
        return False
    return True


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

def cuts_to_observation_vector(cell, cell_cuts, window_radius, n_bins, bin_size=1, take_n_samples=None,
                                   log_distance=False):
    obs = np.zeros(n_bins, dtype=np.int64)

    forward = np.array(list(cell_cuts.keys()))
    if take_n_samples is not None:
        forward = np.random.choice(forward, take_n_samples, replace=True)

    forward.sort()

    total_tests = 0
    for position in forward:
        distance_to_all_points = forward - position
        in_bounds = np.abs(distance_to_all_points[(distance_to_all_points >= -window_radius) & (
                    distance_to_all_points <= window_radius)])
        # Exclude the point itself, which will be of course always associated to a distance 0
        in_bounds = in_bounds[in_bounds > 0] - 1  # Offsets 1bp lower
        total_tests += 1
        # Add 1 to every distance we saw
        if log_distance:
            in_bounds = np.ceil(np.log2(in_bounds) * 100).astype(int)
        else:
            in_bounds = (np.floor(in_bounds / bin_size)).astype(int)
        np.add.at(obs, in_bounds, 1)

    return cell, obs, total_tests


def _cuts_to_observation_vector(kwargs):
        return cuts_to_observation_vector(**kwargs)


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
    argparser.add_argument('-max_distance', type=int,default=2000, help='Maximum distance in both plots and output tables')
    args = argparser.parse_args()

    if not os.path.exists(args.o):
        os.makedirs(args.o)

    # 'Original' analysis
    #analyse(args.alignmentfile, args.o, create_plot=True, verbose=True,strand_specific=False,max_distance=args.max_distance)

    # Stranded analysis:
    sc_cut_dict_stranded = get_sc_cut_dictionary( args.alignmentfile,strand_specific=True,filter_function=strict_read_counts_function)
    distance_counter_fwd_above, distance_counter_fwd_below, distance_counter_rev_above, distance_counter_rev_below = get_stranded_pairwise_counts(sc_cut_dict_stranded)

    # Write tables:
    pd.DataFrame(distance_counter_fwd_above).sort_index().sort_index(1).to_csv(f'{args.o}/STRANDED_fwd_above.csv')
    pd.DataFrame(distance_counter_fwd_below).sort_index().sort_index(1).to_csv(f'{args.o}/STRANDED_fwd_below.csv')
    pd.DataFrame(distance_counter_rev_above).sort_index().sort_index(1).to_csv(f'{args.o}/STRANDED_rev_above.csv')
    pd.DataFrame(distance_counter_rev_below).sort_index().sort_index(1).to_csv(f'{args.o}/STRANDED_rev_below.csv')

    del sc_cut_dict_stranded

    #################
    # Unstranded density analysis:
    sc_cut_dict = get_sc_cut_dictionary( args.alignmentfile,strand_specific=False,filter_function=strict_read_counts_function)
    cpr = get_r1_counts_per_cell(args.alignmentfile)

    def get_commands(one_contig=None):
        for contig in sc_cut_dict:  # sc_cut_dict:
            if '_' in contig or contig in ('chrY', 'chrM', 'chrEBV'):
                continue
            if one_contig is not None and contig != one_contig:
                continue
            for cell, cell_cuts in sc_cut_dict[contig].items():
                yield cell, cell_cuts, contig


    # Calculate distance from one position within a window
    window_radius = args.max_distance
    log_distance = False
    bin_size = 1

    if log_distance:
        n_bins = int(np.ceil(np.log2(window_radius) * 100)) + 1
        x_obs = np.power(2, np.arange(1, n_bins + 1) / 100)

    else:
        n_bins = int(np.ceil(window_radius / bin_size))
        x_obs = np.linspace(1, window_radius + 1, n_bins)  # the associated distance per bin

    # Single cell and one-sided
    # This is a histogram of the amount of observed fragments at distances x:
    obs = defaultdict(lambda: np.zeros(n_bins, dtype=np.int64))
    total_tests = Counter()  # cell -> tests
    with Pool() as workers:
        for cell, cell_obs, n_tests in workers.imap_unordered(
                _cuts_to_observation_vector,

                (
                        {'cell_cuts': cell_cuts,
                         'window_radius': window_radius,
                         'cell': cell,
                         'log_distance': log_distance,
                         'n_bins': n_bins,
                         'bin_size': bin_size,
                         'take_n_samples': None  # sample_target[contig]
                         }
                        for cell, cell_cuts, contig in get_commands()
                )):
            obs[cell] += cell_obs
            total_tests[cell] += n_tests

    p_obs = pd.DataFrame(obs) / pd.Series(total_tests)
    p_obs.index = x_obs

    # Means per library:
    counts = pd.Series(cpr).sort_values()
    window = 35
    df = p_obs.rolling(center=True, window=window).mean()
    df.to_csv(f'{args.o}/strand_unspecific_density.csv')
    df = df[counts[counts > 1_000].index]
    marks = pd.DataFrame({'library': {cell: cell.split('_')[0] for cell in df.columns}})


    fig, ax = plt.subplots(figsize=(15, 8))

    for mark, cells in marks.groupby('library'):
        # if mark!='K9m3':
        #    continue
        df[cells.index].T.iloc[:, 1:].mean(0).iloc[20:].plot(label=f'{mark}, {window}bp window', ax=ax)

    sns.despine()
    ax = plt.gca()
    ax.grid(which='minor')
    ax.grid()
    plt.yscale('log')
    # plt.xscale('log')
    plt.xlabel('distance from cut (bp)')
    plt.ylabel('P(cut)')
    plt.tick_params(axis='y', which='minor')
    # ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    plt.legend()
    plt.savefig(f'{args.o}/density_per_library.png')