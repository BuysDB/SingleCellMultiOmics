import matplotlib
import numpy as np
import pandas as pd
from singlecellmultiomics.utils import is_main_chromosome, get_contig_list_from_fasta
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pysam
import seaborn as sns

# Define chromsome order:
def sort_chromosome_names(l):
    chrom_values = []
    for chrom in l:
        chrom_value = None
        chrom = chrom.replace('chr','').upper()
        if chrom == 'X':
            chrom_value = 99
        elif chrom == 'Y':
            chrom_value = 100
        elif chrom == 'M' or chrom=='MT':
            chrom_value = 101
        elif chrom == 'EBV':
            chrom_value = 102
        elif chrom=='MISC_ALT_CONTIGS_SCMO':
            chrom_value=999
        else:
            chrom_value = int(chrom)
        chrom_values.append(chrom_value)

    indices = sorted(range(len(chrom_values)),key=lambda x:chrom_values[x])
    return [l[idx] for idx in indices]



class GenomicPlot():
    def __init__(self, ref_path, contigs=None):
        """
        Initialise genomic plot

        ref_path(str or pysam.FastaFile) : Path or handle to reference

        """

        if contigs is None:
            self.contigs = sort_chromosome_names(list(filter(is_main_chromosome, get_contig_list_from_fasta(ref_path))))
        else:
            self.contigs = contigs

        # Obtain the lengths:
        if type(ref_path) is str:
            with pysam.FastaFile(ref_path) as reference:
                self.lengths = {r:l for r,l in  zip(reference.references,reference.lengths) if r in self.contigs}
        else:
            self.lengths = {r:l for r,l in  zip(ref_path.references,ref_path.lengths) if r in self.contigs}

        self.total_bp = sum(self.lengths.values())
        # Prune contigs with no length:
        self.contigs = [contig for contig in self.contigs if contig in self.lengths]

    def cn_heatmap(self, df,cell_font_size=3, max_cn=4, method='ward', cmap='bwr',
            figsize=(15,20), xlabel = 'Contigs', ylabel='Cells', **kwargs ):
        """
        Create a heatmap from a copy number matrix

        df: triple indexed dataframe with as columns ('contig', start, end ), as rows cells/samples

        cell_font_size (int): font size of the cell labels

        max_cn (int) : dataframe will be clipped to this value. (Maximum copy number shown)

        method (str) : clustering metric

        cmap (str) : colormap used

        figsize(tuple) : Size of the figure

        xlabel (str) : Label for the x-axis, by default this is Contigs

        ylabel (str) : Label for the x-axis, by default this is Cells

        **kwargs : Arguments which will be passed to seaborn.clustermap

        """

        allelic_mode = len(df.columns[0])==4
        if allelic_mode:
            alleles = [allele for allele in df.columns.get_level_values(0).unique() if not pd.isna(allele)]
            contigs_to_plot = [contig for contig in self.contigs if contig in set(df.columns.get_level_values(1))]
            # Resample the dataframe, drop columns with no allele assigned:
            df = df.loc[:,df.columns.isin(contigs_to_plot, level=1)][alleles].sort_index(1)

            def m(k):
                allele,contig,start,end=k
                return  self.contigs.index(contig), alleles.index(allele),start


            desired_order = sorted( list(df.loc[:,df.columns.isin(self.contigs, level=1)][alleles].sort_index(1).columns), key=m)
            df = df[desired_order]

        else:

            # Figure out what contigs are present in the dataframe:
            contigs_to_plot = [contig for contig in self.contigs if contig in set(df.columns.get_level_values(0))]
            df = df.sort_index(1)[contigs_to_plot]
        try:

            clmap = sns.clustermap(df,
                col_cluster=False,method=method,
                cmap=cmap, vmax=max_cn,vmin=0,
                yticklabels=True, figsize=figsize, **kwargs)
            ax_heatmap = clmap.ax_heatmap
        except Exception as e:
            print(e)
            print('Falling back on heatmap without clustering')

            fig, ax_heatmap = plt.subplots(figsize=figsize)
            clmap = sns.heatmap(df,cmap=cmap,
                vmax=max_cn,vmin=0, yticklabels=True, ax=ax_heatmap, **kwargs)


        prev = None
        xtick_pos = []
        xtick_label = []
        last_idx = 0

        allele = None
        for idx, key in enumerate(df.columns):
            if allelic_mode:
                (allele, contig, start, end)  = key
            else:
                (contig, start, end)  = key

            # Clean up contig label:
            contig = contig.replace('chr', '')
            if allele is not None:
                contig = f'{contig}:{allele}'
            if prev is not None and prev != contig:
                ax_heatmap.axvline(idx-0.5, c='k',lw=1.5, zorder=10)
                xtick_pos.append( (idx+last_idx) / 2)
                xtick_label.append(prev)
                last_idx=idx
            prev = contig

        # Plot last tick..
        xtick_pos.append( (idx+last_idx) / 2)
        xtick_label.append(contig)

        ax_heatmap.set_xticks(xtick_pos)
        ax_heatmap.set_xticklabels(xtick_label,rotation=0, fontsize=8)
        ax_heatmap.set_xlabel(xlabel,labelpad=20)
        ax_heatmap.set_ylabel(ylabel,labelpad=20)

        return clmap

    def get_relative_widths(self):
        return [self.lengths[contig]/self.total_bp for contig in self.contigs]

    def reset_axis(self, contig):
        ax = self[contig]
        ax.clear()
        ax.set_yticklabels([])
        ax.set_yticks([],[])
        ax.set_xticks([],[])
        ax.set_xlabel(contig.replace('chr',''))
        ax.set_xlim(0,self.lengths[contig])

    def get_figure(self, figsize=(20,1)):

        widths = self.get_relative_widths()

        gs_kw = dict(width_ratios=widths)
        figure = plt.figure(figsize =figsize)
        figure.subplots_adjust(bottom=0.25, top=0.75)


        self.gridspec = gridspec.GridSpec(1, len(widths), figure=figure, wspace=0.1, width_ratios=widths)
        self.axis = {}
        prev_ax = None
        for i,contig in enumerate(self.contigs):
           # i = i + 1 # grid spec indexes from 0

            ax = plt.subplot(self.gridspec[i], sharey=prev_ax)
            self.axis[contig] = ax
            self.reset_axis(contig)
            prev_ax=ax
        sns.despine(left=True)
        figure.canvas.draw()
        return figure

    def __getitem__(self, contig):
        return self.axis[contig]
