#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
import singlecellmultiomics

matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
import seaborn as sns
import pysam
import numpy as np
import multiprocessing
from datetime import datetime
from singlecellmultiomics.utils.plotting import GenomicPlot
from singlecellmultiomics.bamProcessing.bamFunctions import verify_and_fix_bam
from singlecellmultiomics.bamProcessing.bamBinCounts import count_fragments_binned, generate_commands, gc_correct_cn_frame, obtain_counts
import os
import argparse
from colorama import Fore, Style


from scipy.cluster.hierarchy import leaves_list
from scipy.cluster.hierarchy import linkage,fcluster
import sklearn.metrics
import random
import collections
import seaborn as sns
import more_itertools
from singlecellmultiomics.utils.pandas import  createRowColorDataFrame
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages

from multiprocessing import Pool
from scipy.optimize import minimize

def square_integer_dist( args,v):
    (amp,offset) = args

    int_distance = np.sum( np.power( (v*amp + offset) - (v*amp + offset).round(),2) )
    distance_from_median =  np.power(np.median((v*amp + offset))-2,2)*30
    return distance_from_median+int_distance

def minimize_square_int_dist(v):
    return minimize(square_integer_dist,
                    x0=[1,0],
                    args=(v, ),
                    bounds=[
                        (0.7,1.3),
                        (-0.01,0.01)]
                   ).x

def minimize_peak_to_integer(df):
    #Calculate amplification factor such that the distance to integers is smallest
    with Pool() as workers:
        minimized = list(workers.imap(minimize_square_int_dist, (row for _,row in df.iterrows() )))

    amps = np.array(minimized)[:,0]
    offsets = np.array(minimized)[:,1]

    return (df.T*amps + offsets).T


def variance_filter(copy_mat, final_segments, segmented_matrix, plot_path=None, vlim=0.025, min_cells_per_seg_call = 5):

    d = copy_mat
    variances = {}
    for chrom, seg in final_segments:
        if not (chrom,seg) in segmented_matrix:
            continue

        median_cn = d[chrom].iloc[:,seg[0]:seg[1]].median(1).round()
        cn_obs = collections.Counter(median_cn).most_common()

        if len(cn_obs)<2:
            continue

        most_common_cn = cn_obs[0][0]
        second_common_cn = cn_obs[1][0]
        if cn_obs[1][1]<min_cells_per_seg_call:
            continue

        variances[chrom,seg] = max( d[chrom].iloc[:,seg[0]:seg[1]].mean(1)[second_common_cn==median_cn].var(0), d[chrom].iloc[:,seg[0]:seg[1]].mean(1)[most_common_cn==median_cn].var(0)) #variance_selected_bins

    vf = pd.DataFrame({'variance':variances})
    variance_selected_bins = vf[ vf['variance']<=vlim ].index

    if plot_path is not None:
        fig, ax = plt.subplots(figsize=(20,4))
        vf.plot.bar(ax=ax)
        ax.axhline(vlim,c='red')
        #vf['variance'].plot.hist()
        #ax = plt.gca()
        #ax.axvline(vlim,c='red')
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.close()

    var_filtered_final_segments = []

    for contig,(start, end) in final_segments:
        k = contig,(start, end)
        if k in variance_selected_bins:
            var_filtered_final_segments.append(k)

    return var_filtered_final_segments


def assign_clusters(copy_mat, final_segments,min_cells_per_cluster=3, MAXCP=4,  min_segment_size = 5):

    d = copy_mat

    segmented_matrix_floating=[]
    segmented_matrix = []
    segmented_matrix_labels = []

    for chrom, seg in final_segments:

        segmented_matrix_floating.append(np.clip(d[chrom].iloc[:,seg[0]:seg[1]].median(1),0, MAXCP  ))
        segmented_matrix.append(np.clip(d[chrom].iloc[:,seg[0]:seg[1]].median(1).round(0),0, MAXCP  ))
        segmented_matrix_labels.append((chrom,seg))

    segmented_matrix = pd.concat(segmented_matrix,1)
    segmented_matrix.columns = segmented_matrix_labels

    segmented_matrix_floating = pd.concat(segmented_matrix_floating,1)
    segmented_matrix_floating.columns = segmented_matrix_labels

    columns_with_info = [len(segmented_matrix[column].unique())>1 for column in segmented_matrix]
    segmented_matrix = segmented_matrix.loc[:,columns_with_info]
    segmented_matrix_labels = np.array(segmented_matrix_labels)[columns_with_info]
    segmented_matrix.columns= pd.MultiIndex.from_tuples(segmented_matrix.columns)

    cnv_clusters = collections.Counter()
    cell_to_unfiltered_cnv = collections.defaultdict(list)
    for cell,row in segmented_matrix.iterrows():
        cnv_clusters[tuple(row)]+=1
        cell_to_unfiltered_cnv[tuple(row)].append(cell)

    keep_clusters = []
    cell_order = []
    cell_cluster_names = []
    current_cluster_name = 1

    median_profiles = []
    print('= Cluster assignment = ')
    total_dropped = 0
    for ci,(cluster, obs) in enumerate( cnv_clusters.most_common() ):
        print(f'\t{cluster}, obs: {obs}')
        if obs>=min_cells_per_cluster:

            keep_clusters.append(cluster)

            cells_in_cluster = []
            for cell in cell_to_unfiltered_cnv[cluster]:
                cells_in_cluster.append(cell)

            # Cluster the cells ..
            cells_in_cluster = np.array(cells_in_cluster)

            cells_in_cluster = cells_in_cluster[leaves_list(
                    linkage(d.loc[cells_in_cluster], method='ward',optimal_ordering=True)
            )]
            #print(obs, len(cells_in_cluster))
            median_profiles.append( d.loc[cells_in_cluster].median(0) )

            cell_order += list(cells_in_cluster)
            cell_cluster_names += [current_cluster_name]*len(cells_in_cluster)

            current_cluster_name+=1
        else:
            print(f'\tDropping cluster {ci}, it has {obs} cells')
            total_dropped+=obs
            #fig, ax = plt.subplots()
            #sns.heatmap(d.loc[cells_in_cluster].sort_index(1)[chrom_order] ,ax=ax, vmax=MAXCP,cmap='bwr')
            #plt.savefig(f'./cluster_{ci}.png')
            #plt.close()


    print(f'{current_cluster_name} clusters identified')
    print(f'{segmented_matrix.shape[0]} cells assigned to a cluster, {total_dropped} cells lost')

    segmented_matrix.columns.names=['contig','range']

    segmented_matrix.index.name='sample'

    return segmented_matrix, cell_cluster_names, median_profiles, cell_order, segmented_matrix_floating


def filter_segment_size(segment_bounds, min_segment_size):
    final_segments = []
    for chrom, bounds_set in segment_bounds.items():
        bounds_list=sorted(list(bounds_set))
        for seg in list( more_itertools.windowed(bounds_list,2) ):
            if np.abs(np.diff(seg))[0] < min_segment_size:
                continue
            seg = (seg[0],seg[1]) # the range is exclusive
            final_segments.append( (chrom,seg) )
    return final_segments

def generate_intitial_clustering(copy_mat, plot_directory, MAXCP=4, chrom_order=None, cn_difference_threshold=0.7, hand_picked_thresholds={}, max_cluster_count=30, seed=42, threshold_offset = 0 ):

    if plot_directory is not None and not os.path.exists(plot_directory):
        os.makedirs(plot_directory)

    if chrom_order is None:
        chrom_order = list(set(copy_mat.columns.get_level_values(0)))

    random.seed(seed)
    segment_bounds = collections.defaultdict(set)
    segment_calls = []
    for chromosome in chrom_order:
        d = copy_mat[chromosome].clip(0,MAXCP)
        L = linkage(d, method='ward')
        scores = []

        if chromosome in hand_picked_thresholds:
            target = hand_picked_thresholds[chromosome]
            print(f'Setting threshold for {chromosome} to {target} clusters')
        else:
            target = 0
        thresholds = list(range(1,max(max_cluster_count,target+1)))
        cluster_count = []
        max_value = None
        max_threshold = None
        max_clustering = None
        clusterings = {}
        for t in thresholds:
            z = fcluster(L, t ,'maxclust')
            clusterings[t] = z
            cluster_count.append(len(set(z)))
            try:
                scores.append( sklearn.metrics.silhouette_score(d,z) )
            except Exception as e:
                scores.append(0)
                pass

            if chromosome in hand_picked_thresholds:
                if t==hand_picked_thresholds[chromosome]:
                    max_value = scores[-1]
                    max_clustering = z
                    max_threshold = t

            elif max_value is None or scores[-1]>max_value :
                max_value = scores[-1]
                max_clustering = z
                max_threshold = t



        print(f'Clustering {chromosome} into {max_threshold} initial clusters')

        max_cluster = max( clusterings.keys() )


        #if max_threshold+1 in clusterings:
    #        max_threshold+=1

        max_threshold+=threshold_offset
        max_threshold=min(max_threshold,max_cluster)

        max_value = scores[max_threshold-1]
        max_clustering = clusterings[max_threshold]


        if plot_directory is not None:
            try:
                plt.plot(thresholds,scores)
                plt.title(chromosome)
                plt.gca().axvline(max_threshold,c='r')
                plt.savefig(f'{plot_directory}/silhouette_score_{chromosome}.png')
                plt.close()
            except Exception as e:
                print(e)

        assignments = max_clustering
        cdf = pd.DataFrame( [assignments],  columns=d.index )
        cdf, lut = createRowColorDataFrame(cdf.T)



        assignments = max_clustering

        delta_cn_hist= collections.Counter()

        if plot_directory is not None:
            sns.clustermap( d.sort_index(1),
                       col_cluster= False, row_cluster=True, method= 'ward', vmax=MAXCP, row_colors=cdf, figsize=(20,40))
            plt.savefig(f'{plot_directory}/clustering_{chromosome}.png')
            plt.close()
            fig, axes = plt.subplots(len(set(assignments)),1,figsize=(8,1 + len(set(assignments))), sharex=True, sharey=True, squeeze=False)
        else:
            axes = [None]*len(set(assignments))

        for ax_col,clust in zip(axes, sorted(list(set(assignments)))):
            if plot_directory is not None:
                ax = ax_col[0]

            data = d[assignments==clust].median().sort_index(0)
            p = 0.005
            sample = data.values
            L = segment(  sample,p=p )

            # Copy number per segment:
            S = validate(sample, L,p=p)

            if plot_directory is not None:
                ax.set_ylim(0, MAXCP+0.5)

            segments = list(more_itertools.windowed(S,2) )

            bps = []
            # @todo: this code has a bug
            for breakpoint, delta_cn in zip(S[1:], np.diff(
                [data.iloc[start:end].median() for start, end in segments] )):
                delta_cn_hist[delta_cn] += 1

                if abs(delta_cn)>cn_difference_threshold:
                    #egment_bounds[chromosome][data.index[min(breakpoint,len(data)-1)]]+=len(data)
                    bps.append(breakpoint)


            bps_including_ends = [1]+bps+[len(data)-2]

            called_segment_indices = list(more_itertools.windowed(bps_including_ends,2) )

            for bp in bps_including_ends:
                #if plot_directory is not None:
                #    ax.axvline( [(start+end)/2 for start, end in data.index.values][bp], c='r')
                segment_bounds[chromosome].add(bp)

            for (start, end) in called_segment_indices:

                data_in_seg = data[start:end]

                color = 'grey'
                if data_in_seg.median()>2.5:
                    color='r'
                if data_in_seg.median()<1.5:
                    color='b'

                if plot_directory is not None:
                    ax.scatter(
                        [(start+end)/2 for start, end in data_in_seg.index.values], data_in_seg.values, s=4,c=color  )



            if plot_directory is not None:
                ax.set_title(f'{chromosome}, cluster {clust}')
                ax.set_ylabel('copy number')
                sns.despine(ax=ax)
                for cp in range(1,MAXCP+1):
                    ax.axhline(cp,c='k',lw=0.5)


        if plot_directory is not None:
            ax.set_xlabel('position')
            plt.tight_layout()

            plt.savefig(f'{plot_directory}/segments_{chromosome}.png',dpi=120)
            plt.close()

    return segment_bounds





def cbs_stat(x):
    '''Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t.
    Returns t, i0, i1'''
    x0 = x - np.mean(x)
    n = len(x0)
    y = np.cumsum(x0)
    e0, e1 = np.argmin(y), np.argmax(y)
    i0, i1 = min(e0, e1), max(e0, e1)
    s0, s1 = y[i0], y[i1]
    return (s1-s0)**2*n/(i1-i0+1)/(n+1-i1+i0), i0, i1+1


def tstat(x, i):
    '''Return the segmentation statistic t testing if i is a (one-sided)  breakpoint in x'''
    n = len(x)
    s0 = np.mean(x[:i])
    s1 = np.mean(x[i:])
    return (n-i)*i/n*(s0-s1)**2






def cbs(x, shuffles=1000, p=.05):
    '''Given x, find the interval x[i0:i1] with maximal segmentation statistic t. Test that statistic against
    given (shuffles) number of random permutations with significance p.  Return True/False, t, i0, i1; True if
    interval is significant, false otherwise.'''

    max_t, max_start, max_end = cbs_stat(x)
    if max_end-max_start == len(x):
        return False, max_t, max_start, max_end
    if max_start < 5:
        max_start = 0
    if len(x)-max_end < 5:
        max_end = len(x)
    thresh_count = 0
    alpha = shuffles*p
    xt = x.copy()
    for i in range(shuffles):
        np.random.shuffle(xt)
        threshold, s0, e0 = cbs_stat(xt)
        if threshold >= max_t:
            thresh_count += 1
        if thresh_count > alpha:
            return False, max_t, max_start, max_end
    return True, max_t, max_start, max_end


def rsegment(x, start, end, L=[], shuffles=1000, p=.05):
    '''Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
    '''
    threshold, t, s, e = cbs(x[start:end], shuffles=shuffles, p=p)
    if (not threshold) | (e-s < 5) | (e-s == end-start):
        L.append((start, end))
    else:
        if s > 0:
            rsegment(x, start, start+s, L)
        if e-s > 0:
            rsegment(x, start+s, start+e, L)
        if start+e < end:
            rsegment(x, start+e, end, L)
    return L


def segment(x, shuffles=1000, p=.05):
    '''Segment the array x, using significance test based on shuffles rearrangements and significance level p
    '''
    start = 0
    end = len(x)
    L = []
    rsegment(x, start, end, L, shuffles=shuffles, p=p)
    return L


def validate(x, L, shuffles=1000, p=.01):
    S = [x[0] for x in L]+[len(x)]
    SV = [0]
    left = 0
    for test, s in enumerate(S[1:-1]):
        t = tstat(x[S[left]:S[test+2]], S[test+1]-S[left])
        threshold = 0
        thresh_count = 0
        site = S[test+1]-S[left]
        xt = x[S[left]:S[test+2]].copy()
        flag = True
        for k in range(shuffles):
            np.random.shuffle(xt)
            threshold = tstat(xt, site)
            if threshold > t:
                thresh_count += 1
            if thresh_count >= p*shuffles:
                flag = False
                break
        if flag:
            SV.append(S[test+1])
            left += 1
    SV.append(S[-1])
    return SV


def get_segment_calls(copy_series, p_value=0.01,  shuffles=10000, plot=None):

    calls = {}
    for contig in set(copy_series.columns.get_level_values(0)):
        calls[contig] = []
        copy_vector = copy_series[contig].median().sort_index()

        x = [ 0.5*(s+e) for s,e in copy_vector.index ]

        if plot is not None:
            plot.axis[contig].scatter(x, copy_vector.values, c='grey', s=1)

        l = segment( copy_vector.values, p=p_value, shuffles=shuffles)

        s_calls = validate(copy_vector.values, l, p=p_value, shuffles=shuffles)

        if s_calls[-1]!=len(copy_vector.values):
            s_calls.append(len(copy_vector.values))

        segments = list(more_itertools.windowed(s_calls,2) )

        prev_cn = None
        for segidx, (start,end) in enumerate(segments):
            y = copy_vector.values[start:end]

            try:
                cn = int(np.round(np.median(y)))
            except ValueError:
                cn = 0 if prev_cn is None else prev_cn


            if prev_cn is not None and prev_cn==cn:
                continue

            calls[contig].append( (start, end) )

            # Obtain the actual coordinates of the segment..
            #print(y,)

            x = []
            for i, (s,e) in enumerate(copy_vector.index[start:end]):

                x.append((s + e) * 0.5)
                """
                if i==0:
                    x.append(s)
                else:
                    if i==(end-start-1):
                        x.append( e )
                    else:
                        x.append( (s+e)*0.5 )
                """

            if plot is not None and contig in plot.axis:
                color = {0:'k',1:'b',2:'grey',3:'red',4:'orange',5:'cyan'}.get(cn,'cyan')
                plot.axis[contig].scatter(x, y, c=color, s=2)

            if prev_cn is not None and plot is not None and contig in plot.axis:

                #plot.axis[contig].axvline(copy_vector.index[start][0],c='k', lw=1)
                if False:
                    try:
                        plot.axis[contig].axvline(copy_vector.index[end][1],c='k', lw=1)
                    except IndexError:
                        plot.axis[contig].axvline(copy_vector.index[end-1][1],c='k', lw=1)

            prev_cn=cn
            #plot.axis[contig].axvline(copy_vector.index[end-1],c='k')

    return calls

def bulk_trace(pdf_path, copy_mat, cell_cluster_names, cell_order,segmented_matrix_floating,segmented_matrix ):
    with PdfPages(pdf_path,
                  metadata={'Creator': f'SingleCellMultiOmics {singlecellmultiomics.__version__}', 'Author': 'SCMO',
                            'Title': 'Traces for assigned clusters'}) as pdf:


        for cluster in set(cell_cluster_names):
            fig = h.get_figure()
            cells_in_cluster = np.array(cell_order)[np.array(cell_cluster_names)==cluster]
            get_segment_calls(copy_mat.loc[cells_in_cluster][chrom_order], plot=h, p_value=0.01,  shuffles=10000)
            plt.suptitle(f'Cluster {cluster}, {len(cells_in_cluster)} cells,  {len(cells_in_cluster)*100/len(cell_cluster_names):.2f}% of total')
            pdf.savefig(fig)
            plt.close()

        # Create histograms per segment
        for segment in segmented_matrix:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.hist(segmented_matrix_floating[segment], bins=40)
            sns.despine()
            ax.set_xlabel('copy number')
            ax.set_ylabel('# cells')
            plt.title(f'segment {segment[0]}:{segment[1][0] * bin_size}-{segment[1][1] * bin_size}', pad=30)
            ax.set_xlim(0, MAXCP + 0.5)
            for boundary in np.arange(0.5, MAXCP, 1):
                plt.axvline(boundary, c='k', lw=1)

            # Write total cells in CN:
            # ax.set_ylim(0,ax.get_ylim()[1]*1.1)
            y = ax.get_ylim()[1]
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            for cn in range(0, MAXCP + 1):
                n_cells = (segmented_matrix[segment] == cn).sum()
                ax.text(cn + (0.25 if cn == 0 else 0), y, f'CN:{cn}\n{n_cells} cells', horizontalalignment='center')
            ax.grid(axis='y', which='both')
            pdf.savefig(fig)
            plt.close()


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Export and plot copy number profiles
    """)
    argparser.add_argument('bamfiles', metavar='bamfiles', type=str, nargs='+')
    argparser.add_argument('-ref', help='path to reference fasta', type=str, required=True)
    argparser.add_argument('-bin_size', default=500_000, type=int)
    argparser.add_argument('-max_cp', default=5, type=int)
    argparser.add_argument('-threads', default=16, type=int)
    argparser.add_argument('-bins_per_job', default=5, type=int)
    argparser.add_argument('-pct_clip', default=99.999, type=float)
    argparser.add_argument('-min_mapping_qual', default=40, type=int)
    argparser.add_argument('-molecule_threshold', default=5_000, type=int)
    argparser.add_argument('-ignore_contigs', default=None, type=str, help='Comma separated contigs to ignore for the analysis')




    argparser.add_argument('--ignore_mp',action='store_true',help='Ignore mp tag value')
    argparser.add_argument('--ignore_qcfail',action='store_true',help='Ignore qcfail tag value')
    argparser.add_argument('--allelic',action='store_true',help='Perform allele specific analysis (requires DA tag)')
    argparser.add_argument('-rawmatplot', type=str, help='Path to raw matrix, plot is not made when this path is not supplied ')
    argparser.add_argument('-gcmatplot', type=str, help='Path to gc corrected matrix, plot is not made when this path is not supplied ')
    argparser.add_argument('-histplot', type=str, help='Path to histogram ')

    argparser.add_argument('-rawmat', type=str)
    argparser.add_argument('-countmat', type=str)
    argparser.add_argument('-gcmat', type=str)

    argparser.add_argument('-norm_method', default='median', type=str, help='Either mean or median')


    cops = argparser.add_argument_group('Clustering options')

    cops.add_argument('-clustering_output_folder', type=str)
    cops.add_argument('-min_segment_size', default=5, type=int)
    cops.add_argument('-min_cells_per_cluster', default=8, type=int)
    cops.add_argument('-vlim', default=0.04, type=float, help= 'variance limit')
    cops.add_argument('-cn_difference_threshold', default=0.7, type=float)
    cops.add_argument('-n_clusters_for_contig', default=None, type=str,
        help="""Use to manually set the amount of clusters present at a contig
            when the algorithm is underclustering.
            Format: chr1:10,chr3:7
            , when not supplied the silhouette_score is used to determine the amount of clusters""")

    args = argparser.parse_args()

    alignments_paths = args.bamfiles
    bin_size = args.bin_size
    MAXCP = args.max_cp
    pct_clip = args.pct_clip
    bins_per_job = args.bins_per_job
    min_mapping_qual = args.min_mapping_qual
    threads = args.threads
    molecule_threshold = args.molecule_threshold

    histplot=args.histplot
    rawmatplot=args.rawmatplot
    gcmatplot=args.gcmatplot
    rawmat=args.rawmat
    gcmat=args.gcmat

    ignore_contigs = None if args.ignore_contigs is None else args.ignore_contigs.split(',')

    reference = pysam.FastaFile(args.ref)
    h=GenomicPlot(reference, ignore_contigs=ignore_contigs)
    contigs = GenomicPlot(reference).contigs

    kwargs = {'ignore_mp':args.ignore_mp,'ignore_qcfail':args.ignore_qcfail}
    print("Creating count matrix ... ")

    # Check if the bam files are in good shape:
    for path in alignments_paths:
        verify_and_fix_bam(path)

    if args.allelic:
        key_tags=['DA']
    else:
        key_tags=None

    commands = generate_commands(
                alignments_paths,
                bin_size=bin_size,key_tags=key_tags,
                bins_per_job=bins_per_job,head=None,min_mq=min_mapping_qual,kwargs=kwargs )

    counts = obtain_counts(commands,
                            reference=reference,
                            threads=threads,
                            live_update=False,
                            show_n_cells=None,
                            update_interval=None,show_progress=True )

    print(f"Creating count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")



    if histplot is not None:
        print("Creating molecule histogram ... ",end="")
        df = pd.DataFrame(counts).T.fillna(0)
        fig, ax = plt.subplots()
        cell_sums = df.sum()
        cell_sums.name = 'Frequency'
        cell_sums.plot.hist(bins=50)
        ax.set_xlabel('# molecules')
        ax.set_xscale('log')
        ax.axvline(molecule_threshold, c='r', label='Threshold')
        plt.legend()
        plt.savefig(histplot)
        print(f"\rCreating molecule histogram [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")
        plt.close()


    # Convert the count dictionary to a dataframe

    df = pd.DataFrame(counts).T.fillna(0)

    if df.shape[0]==0:
        raise ValueError('Resulting count matrix is empty. Is this file correctly tagged? Try adding the --ignore_mp flag')

    if args.countmat is not None:
        print("Exporting count matrix ... ", end="")
        if args.countmat.endswith('.pickle.gz'):
            df.to_pickle(args.countmat)
        else:
            df.to_csv(args.countmat)

        print(f"\rExporting count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

    if args.allelic:
        alleles = [allele for allele in df.index.get_level_values(0).unique() if not pd.isna(allele)]
        print(f'Found alleles: {", ".join(alleles)}')
        df = df.loc[alleles]

    print("Filtering count matrix ... ", end="")
    # remove cells were the median is zero
    if args.norm_method=='median':
        try:
            shape_before_median_filter = df.shape
            df = df.T[df.median()>0].T
            shape_after_median_filter = df.shape
            print(shape_before_median_filter,shape_after_median_filter )
            # Remove rows with little counts
            df = df.T[df.sum()>molecule_threshold].T
            df = df / np.percentile(df,pct_clip,axis=0)
            df = np.clip(0,MAXCP,(df / df.median())*(2 if not args.allelic else 1))
            df = df.T
        except Exception as e:
            print(f"\rMedian normalisation [ {Fore.RED}FAIL{Style.RESET_ALL} ] ")
            args.norm_method = 'mean'

    if args.norm_method == 'mean':
        shape_before_median_filter = df.shape
        df = df.T[df.mean()>0].T
        shape_after_median_filter = df.shape
        # Remove rows with little counts
        df = df.T[df.sum()>molecule_threshold].T
        df = df / np.percentile(df,pct_clip,axis=0)
        df = np.clip(0,MAXCP,(df / df.mean())* (2 if not args.allelic else 1))
        df = df.T

    if args.norm_method == 'spikein':
        print(df)
        df =  (df.T/(df['J02459.1'].sum()) ).T

        df = np.clip(0,MAXCP,df)

    # Perform peak transfromation to ensure integer copy numbers with median of ~2 :
    if args.norm_method=='median':
        try:
            df = minimize_peak_to_integer(df)
        except Exception as e:
            print("minimize_peak_to_integer, failed:")
            print(e)

    if df.shape[0]==0:
        print(f"\rRaw count matrix [ {Fore.RED}FAIL{Style.RESET_ALL} ] ")
        raise ValueError('Resulting count matrix is empty, review the filter settings')
    else:
        print(f"\rFiltering count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")
        print( f'{df.shape[0]} cells, and {df.shape[1]} bins remaining' )
    del counts

    if rawmat is not None:
        print("Exporting raw count matrix ... ", end="")
        if rawmat.endswith('.pickle.gz'):
            df.to_pickle(rawmat)
        else:
            df.to_csv(rawmat)

        print(f"\rExporting raw count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

    if rawmatplot is not None:
        print("Creating raw heatmap ...", end="")
        h.cn_heatmap(df, figsize=(15*(2 if args.allelic else 1),15+(0.05)*df.shape[0]))
        plt.savefig(rawmatplot)
        print(f"\rCreating raw heatmap [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")
        plt.close('all')

    if gcmatplot is not None or gcmat is not None or args.clustering_output_folder is not None:
        print("Performing GC correction ...", end="")
        corrected_cells = gc_correct_cn_frame(df, reference, MAXCP, threads, norm_method=args.norm_method)
        print(f"\rPerforming GC correction [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

        # Perform peak transform to ensure integer copy numbers with median of ~2 :
        if args.norm_method=='median':
            try:
                corrected_cells = minimize_peak_to_integer(corrected_cells)
            except Exception as e:
                print("minimize_peak_to_integer, failed:")
                print(e)


    if gcmatplot is not None:
        print("Creating heatmap ...", end="")
        h.cn_heatmap(corrected_cells,figsize=(15*(2 if args.allelic else 1),15))
        plt.savefig(gcmatplot)
        plt.close('all')
        print(f"\rCreating heatmap [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")

    if gcmat is not None:
        print("Exporting corrected count matrix ... ")
        if gcmat.endswith('.pickle.gz'):
            corrected_cells.to_pickle(gcmat)
        else:
            corrected_cells.to_csv(gcmat)
        print(f"\rExporting corrected count matrix [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")


    if args.clustering_output_folder is not None:

        clustering_plot_folder = f'{args.clustering_output_folder}/plots'
        if not os.path.exists( clustering_plot_folder ):
            os.makedirs(clustering_plot_folder)

        segmentation_plot_folder = f'{args.clustering_output_folder}/plots/initial_segmentation'
        if not os.path.exists( segmentation_plot_folder ):
            os.makedirs(segmentation_plot_folder)


        copy_mat = corrected_cells
        chrom_order = [c for c in h.contigs]
        min_cells_per_cluster = args.min_cells_per_cluster
        min_segment_size = args.min_segment_size # segment size in bins
        hand_picked_thresholds = {}

        print("Creating first rough clustering")
        if args.n_clusters_for_contig is not None:


            hand_picked_thresholds = {
                contig_n.split(':')[0]:int(contig_n.split(':')[1])
                for contig_n in args.n_clusters_for_contig.split(',')
                }

            for c in hand_picked_thresholds:
                if c not in chrom_order:
                    print(f"The contig {c} of which a threshold of {hand_picked_thresholds[c]} was set, is not part of the clustering. Make sure the name of the contig is correct. For example, it could be that the supplied reference does not use a chr prefix. Pick from: {chrom_order[:30]} ..." )

        segment_bounds = generate_intitial_clustering(copy_mat,
                             plot_directory=segmentation_plot_folder,
                             MAXCP=MAXCP,
                             chrom_order=chrom_order,
                             hand_picked_thresholds=hand_picked_thresholds,
                             cn_difference_threshold=args.cn_difference_threshold,


                             seed=42 )

        # Filter for segment size
        print("Filtering for segment size")
        segments = filter_segment_size(segment_bounds, min_segment_size=min_segment_size )
        segmented_matrix, cell_cluster_names, median_profiles,cell_order, segmented_matrix_floating = assign_clusters(copy_mat,
                                                                                  segments,
                                                                                  MAXCP=MAXCP,
                                                                                   min_cells_per_cluster=min_cells_per_cluster,
                                                                                   min_segment_size=min_segment_size)
        bulk_trace(f'{clustering_plot_folder}/segments_wo_variance_filter.pdf', copy_mat, cell_cluster_names, cell_order,segmented_matrix_floating, segmented_matrix)

        cell_annot_df = pd.DataFrame([cell_cluster_names, [cell.split('_')[0] for cell in cell_order]],
                                     columns=cell_order).T
        cell_annot_colors, lut = createRowColorDataFrame(cell_annot_df)
        cell_annot_colors.columns=['cluster','library']
        cell_annot_df.columns=['cluster','library']
        h.cn_heatmap(copy_mat.sort_index(1)[chrom_order].loc[cell_order],
                     row_colors=cell_annot_colors,
                     figsize=(30,30), row_cluster=False)
        plt.savefig(f'{clustering_plot_folder}/segment_based_clustering_wo_variance_filter.png',dpi=150)
        cell_annot_df.to_csv(f'{args.clustering_output_folder}/cell_clusters_wo_variance_filter.csv')
        cell_annot_df.to_pickle(f'{args.clustering_output_folder}/cell_clusters_wo_variance_filter.pickle.gz')
        segmented_matrix.to_csv(f'{args.clustering_output_folder}/segmented_matrix_wo_variance_filter.csv')
        segmented_matrix.to_pickle(f'{args.clustering_output_folder}/segmented_matrix_wo_variance_filter.pickle.gz')

        # Variance filter
        print("Variance filter")
        var_filtered_final_segments = variance_filter(copy_mat, segments, segmented_matrix, f'{clustering_plot_folder}/segment_variance.png',vlim=args.vlim)

        print("Creating final segmentation")
        segmented_matrix_f, cell_cluster_names, median_profiles, cell_order, segmented_matrix_floating = assign_clusters(copy_mat,
                                                                                              var_filtered_final_segments,
                                                                                              min_cells_per_cluster=min_cells_per_cluster,
                                                                                              MAXCP=MAXCP, min_segment_size=min_segment_size)
        # Create plot of clustering:
        print("Creating plots and tables")
        cell_annot_df = pd.DataFrame([cell_cluster_names, [cell.split('_')[0] for cell in cell_order]],
                                     columns=cell_order)
        cell_annot_df.to_csv(f'{args.clustering_output_folder}/cell_clusters.csv')
        cell_annot_df.to_pickle(f'{args.clustering_output_folder}/cell_clusters.pickle.gz')
        cell_annot_df, lut = createRowColorDataFrame(cell_annot_df.T)
        cell_annot_df.columns=['cluster','library']

        h.cn_heatmap(copy_mat.sort_index(1)[chrom_order].loc[cell_order],
                     row_colors=cell_annot_df,
                     figsize=(30,30), row_cluster=False)
        plt.savefig(f'{clustering_plot_folder}/segment_based_clustering.png',dpi=150)
        plt.close()

        sns.clustermap( segmented_matrix_f.loc[cell_order], row_colors=cell_annot_df, row_cluster=False, col_cluster=False )
        plt.savefig(f'{clustering_plot_folder}/segment_clustering.png',dpi=150)
        plt.close()
        segmented_matrix_f.to_csv(f'{args.clustering_output_folder}/segmented_matrix.csv')
        segmented_matrix_f.to_pickle(f'{args.clustering_output_folder}/segmented_matrix.pickle.gz')

        # Create bulk trace plot:
        bulk_trace(f'{clustering_plot_folder}/segments.pdf', copy_mat, cell_cluster_names, cell_order,segmented_matrix_floating,segmented_matrix_f)
