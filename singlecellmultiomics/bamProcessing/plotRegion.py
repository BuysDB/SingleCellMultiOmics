#!/usr/bin/env python


from multiprocessing import Pool
from singlecellmultiomics.bamProcessing.bamFunctions import mate_iter
import argparse
import pysam
from glob import glob
import pandas as pd
from singlecellmultiomics.bamProcessing import get_contig_sizes
from collections import Counter, defaultdict
from singlecellmultiomics.features import FeatureContainer
import os
from matplotlib.patches import Rectangle
import matplotlib as mpl
from scipy.ndimage import gaussian_filter
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from singlecellmultiomics.bamProcessing import get_contigs_with_reads





def _generate_count_dict(args):

    bam_path, bin_size, contig, start, stop = args #reference_path  = args

    #reference_handle = pysam.FastaFile(reference_path)
    #reference = CachedFasta(reference_handle)


    cut_counts = defaultdict(Counter )
    i = 0
    with pysam.AlignmentFile(bam_path) as alignments:

        for R1,R2 in mate_iter(alignments, contig=contig, start=start, stop=stop):

            if R1 is None or R1.is_duplicate or not R1.has_tag('DS') or R1.is_qcfail:
                continue

            cut_pos = R1.get_tag('DS')
            sample = R1.get_tag('SM')

            bin_idx=int(cut_pos/bin_size)*bin_size
            cut_counts[(contig,bin_idx)][sample] += 1

    return cut_counts, contig, bam_path


def get_binned_counts(bams, bin_size, regions=None):

    fs = 1000
    if regions is None:
        regions = [(c,None,None) for c in get_contig_sizes(bams[0]).keys()]

    else:
        for i,r in enumerate(regions):
            if type(r)==str:
                regions[i] = (r,None,None)
            else:
                contig, start, end =r
                if type(start)==int:
                    start = max(0,start-fs)

                regions[i] = (contig,start,end)

    jobs = [(bam_path, bin_size, *region) for region, bam_path in product(regions, bams)]


    cut_counts = defaultdict(Counter)
    with Pool() as workers:

        for i, (cc, contig, bam_path) in enumerate(workers.imap(_generate_count_dict,jobs)):

            for k,v in cc.items():
                cut_counts[k] += v

            print(i,'/', len(jobs), end='\r')

    return pd.DataFrame(cut_counts).T



def plot_region(counts, features, contig, start, end, sigma=2, target=None, caxlabel='Molecules per spike-in'):

    if target is None:
        target = f'{contig}_{start}_{end}.png'

    def create_gene_models(start,end,ax):

        exon_height = 0.010
        gene_height = 0.0002
        spacer = 0.035

        overlap_dist = 200_000


        gene_y = {}
        ymax = 0
        for fs,fe,name,strand, feature_meta in features.features[contig]:


            if not (((fs>=start or fe>=start) and (fs<=end or fe<=end))):
                continue
            feature_meta = dict(feature_meta)

            if feature_meta.get('type') == 'gene':

                if not 'gene_name' in  feature_meta or feature_meta.get('gene_name').startswith('AC'):
                    continue

                # Determine g-y coordinate:

                gy_not_avail = set()
                for gene,(s,e,loc) in gene_y.items():
                    if (s+overlap_dist>=fs and s-overlap_dist<=fe) or (e+overlap_dist>=fs and e-overlap_dist<=fe):
                        # Overlap:
                        gy_not_avail.add(loc)

                gy = 0
                while gy in gy_not_avail:
                    gy+=1


                gene_y[name] = (fs,fe,gy)

                y_offset = gy * spacer

                ymax = max(y_offset+gene_height,ymax)

                r = Rectangle((fs,-gene_height*0.5 + y_offset), fe-fs, gene_height, angle=0.0, color='k')

                ax.add_patch( r )
                ax.text((fe+fs)*0.5,-1.6*exon_height + y_offset,feature_meta.get('gene_name'),horizontalalignment='center',
              verticalalignment='center',fontsize=3)
                #print(feature_meta)



        if False:

            for xx in range(3):
                for fs,fe,name,strand, feature_meta in features.features[contig]:

                    if not (((fs>=start or fe>=start) and (fs<=end or fe<=end))):
                        continue

                    feature_meta = dict(feature_meta)
                    if not name in gene_y:
                        continue

                    if feature_meta.get('type') == 'exon':
                        y_offset = gene_y[name][2]*spacer
                        ymax = max(y_offset+exon_height,ymax)
                        r = Rectangle((fs,-exon_height*0.5 + y_offset), fe-fs, exon_height, angle=0.0,color='k', lw=0)
                        ax.add_patch( r )



        ax.set_xlim(start,end)
        ax.set_ylim(-0.1,ymax)
        #ax.axis('off')
        ax.set_yticks([])
        ax.set_xlabel(f'chr{contig} location bp', fontsize=6)

        #print([t.get_text() for t in ax.get_xticklabels()])
        #ax.set_xticklabels([t.get_text() for t in ax.get_xticklabels()],fontsize=4)
        ax.set_xticklabels(ax.get_xticks(), fontsize=4)


        ax.tick_params(length=0.5)


    for sigma in range(2,3):

        mpl.rcParams['figure.dpi'] = 300

        font = {'family' : 'helvetica',
                'weight' : 'normal',
                'size'   : 8}

        mpl.rc('font', **font)

        if end - start < 3_000_000:
            mode ='k'
            stepper = 100_000
            res = 100
        else:
            mode='M'
            stepper=1_000_000
            res = 1



        qf = counts.loc[:, [(c,p) for c,p in counts if c==contig and p>=start and p<=end] ].sort_index()
        qf = qf.sort_index(1).sort_index(0)
        qf = pd.DataFrame(gaussian_filter(qf, sigma=(0.00001,sigma)), index=qf.index, columns=qf.columns)
        qf = qf.sort_index(1).sort_index(0)

        cm = sns.clustermap(qf,
               #z_score=0,
                row_cluster=False,
                col_cluster=False,
                vmax=np.percentile(qf,99.5),#0.0005,
                #vmax=10,
                dendrogram_ratio=0.1,
                #row_colors=row_colors.loc[qf.index].drop('LOWESS_STAGE',1),
                figsize=(8,4), cmap='Greys', cbar_kws={"shrink": .1},
                cbar_pos=(0.0, 0.5, 0.01, 0.16),)

        ax = cm.ax_col_dendrogram
        qf.mean().plot.bar(ax=ax,color='k',width=1)


        ax.set_yticks([])

        cm.ax_heatmap.set_xticks([]) #np.arange(start,end, 1_000_000))
        cm.ax_heatmap.set_yticks([])


        cm.ax_heatmap.set_ylabel(f'{qf.shape[0]} single cells', fontsize=8)
        cm.ax_heatmap.tick_params(length=0.5)
        cm.ax_heatmap.set_xlabel(None)

        ax.grid()
        cm.cax.set_ylabel(caxlabel,fontsize=6)
        cm.cax.tick_params(labelsize=4)
        #plt.suptitle(mark, x=0.05)


        fig = plt.gcf()
        heatmap_start_x,heatmap_start_y, heatmap_end_x, heatmap_end_y = cm.ax_heatmap.get_position().bounds

        width = heatmap_end_x #-heatmap_start_x
        height = 0.2 if features is not None else 0.05
        ax = fig.add_axes(  (heatmap_start_x, heatmap_start_y-height-0.02, width, height)  )
        ax.ticklabel_format(axis='x',style='sci')

        sns.despine(fig=fig, ax=ax)
        if features is not None:
            create_gene_models(start,end,ax=ax)
        else:
            ax.set_xlim(start,end)
            #ax.axis('off')
            ax.set_yticks([])
            ax.set_xlabel(f'chr{contig} location bp', fontsize=6)

            #ax.set_xticklabels(ax.get_xticks(), fontsize=4)
            plt.xticks(fontsize=4)
            ax.tick_params(length=0.5)

        plt.savefig(target)
        plt.close()


if __name__=='__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Plot a genomic region')

    argparser.add_argument('bams', type=str, nargs='+', help='(X) Training bam files')

    argparser.add_argument('-regions', type=str, help='Regions to plot, with a bin size behind it, for example: 1:1000-100000:1000 , will be a single region plotted with a 1000bp bin size split regions by commas without a space')
    argparser.add_argument('-features', type=str, help='Gene models to plot (.gtf file or .gtf.gz)', required=False)
    argparser.add_argument('-norm', type=str, help='Normalize to, select from : total-molecules,spike', default='total-molecules')
    argparser.add_argument('-prefix', type=str, help='Prefix for output file',default='')
    argparser.add_argument('-format', type=str, help='png or svg',default='png')

    args = argparser.parse_args()

    regions = []
    contigs = set()
    for region in args.regions.split(','):
        contig = region.split(':')[0]
        if not '-' in region:
            start, end = None, None
        else:
            start, end = region.split(':')[1].split('-')
            start = int(start)
            end = int(end)
        bin_size = int(region.split(':')[-1])
        if start is not None:
            print(f'Region: {contig} from {start} to {end} with bin size : {bin_size}')
        else:
            print(f'Region: {contig} with bin size : {bin_size}')
        contigs.add(contig)
        regions.append( ((contig,start,end), bin_size))

    contigs=list(contigs)
    bams = args.bams

    if args.features is not None:
        print('Reading features')
        features = FeatureContainer()
        if len(contigs)==1:
            print(f'Reading only features from {contigs[0]}')
            features.loadGTF(args.features,store_all=True,contig=contigs[0])
        else:
            features.loadGTF(args.features,store_all=True)
    else:
        features = None
    print('Counting')

    # Obtain counts per cell
    norm = 'spike'
    if norm == 'spike':
        normalize_to_counts = get_binned_counts(bams, bin_size=10_000_000, regions=['J02459.1'])
    elif norm=='total-molecules':
        normalize_to_counts = get_binned_counts(bams, bin_size=10_000_000)

    for region, region_bin_size in regions:
        print(f'Plotting {region}')
        contig, start, end = region
        region_counts = get_binned_counts(bams, region_bin_size, regions=[ region ] )
        counts = (region_counts/normalize_to_counts.sum()).fillna(0).T.sort_index(1).sort_index(0)

        # Fill non intialized bins with zeros:
        add = []
        for i in np.arange(counts.columns[0][1], counts.columns[-1][1], region_bin_size):
            if not (contig,i) in counts.columns:
                add.append((contig,i))

        for a in add:
            counts[a] = 0

        counts = counts.sort_index(1)


        target = args.prefix+f'{contig}_{start}-{end}_{region_bin_size}.{args.format}'
        plot_region(counts, features, contig, start, end, sigma=2, target=target, caxlabel='Molecules per spike-in' if norm =='spike' else 'Molecules / total molecules')
