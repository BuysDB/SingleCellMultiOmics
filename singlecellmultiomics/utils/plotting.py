import matplotlib
import numpy as np
import pandas as pd
from singlecellmultiomics.utils import is_main_chromosome, get_contig_list_from_fasta
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pysam
import seaborn as sns
from matplotlib.patches import Circle
from itertools import product
import collections
import string
import math
import matplotlib.patches as mpatches


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
            try:
                chrom_value = int(chrom)
            except Exception as e:
                chrom_value = 999 + sum((ord(x) for x in chrom))
        chrom_values.append(chrom_value)

    indices = sorted(range(len(chrom_values)),key=lambda x:chrom_values[x])
    return [l[idx] for idx in indices]



class GenomicPlot():
    def __init__(self, ref_path, contigs=None, ignore_contigs=None):
        """
        Initialise genomic plot

        ref_path(str or pysam.FastaFile) : Path or handle to reference

        """

        if contigs is None:
            self.contigs = sort_chromosome_names(list(filter(lambda x: is_main_chromosome(x) and (ignore_contigs is None or x not in ignore_contigs) , get_contig_list_from_fasta(ref_path))))
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

    def cn_heatmap(self, df,cell_font_size=3, max_cn=4, method='ward', cmap='bwr', yticklabels=True,
            figsize=(15,20), xlabel = 'Contigs', ylabel='Cells', vmin=0, xtickfontsize=8, **kwargs ):
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
                yticklabels=yticklabels, figsize=figsize, **kwargs)
            ax_heatmap = clmap.ax_heatmap
        except Exception as e:
            print(e)
            print('Falling back on heatmap without clustering')

            fig, ax_heatmap = plt.subplots(figsize=figsize)
            clmap = sns.heatmap(df,cmap=cmap,
                vmax=max_cn,vmin=vmin, yticklabels=True, ax=ax_heatmap, **kwargs)


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
        ax_heatmap.set_xticklabels(xtick_label,rotation=0, fontsize=xtickfontsize)
        ax_heatmap.set_xlabel(xlabel,labelpad=20)
        ax_heatmap.set_ylabel(ylabel,labelpad=20)

        return clmap

    def get_relative_widths(self):
        return [self.lengths[contig]/self.total_bp for contig in self.contigs]

    def reset_axis(self, contig):
        ax = self[contig]
        ax.clear()
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_xticks([])
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


def plot_plate_layout(plate_layout:dict , welllabel2coord:dict , suptitle="Plate layout", cmap=None,class_colors=None) -> dict:
    """
    Plate layout is a dictionary of:
    'D10': 'condition A'
    'D11': 'condition B',
    etc.
    use "empty" for empty wells

    welllabel2coord:
    {'A0': (0, 0),
     'A1': (0, 1),
     'A2': (0, 2),
     'A3': (0, 3), .. }
    """
    states = list(set(plate_layout.values()).difference(set(['empty'])))
    state_to_index = {state: i for i,state in enumerate(states)}
    state_to_index['empty'] = np.nan
    plate_indices = {welllabel2coord[well]:state_to_index[value] for well,value in plate_layout.items()}

    if class_colors is not None:
        plate_values = {welllabel2coord[well]:value for well,value in plate_layout.items()}
        fig,ax,cbar = plot_plate(plate_values,class_colors=class_colors)
        return {'fig':fig,'ax':ax,'cbar':cbar,'state_to_index':state_to_index, 'plate_indices':plate_indices}

    if cmap is None:
        cmap = plt.get_cmap('tab10').copy()
        cmap.set_bad( (0.2,0.2,0.2) )
        cmap.set_under(color='k')


    fig,ax,cbar = plot_plate(plate_indices,vmax=cmap.N,cmap=cmap,log=False,vmin=0,usenorm=True,colorbarargs={'extend':'min'})

    cbar.set_yticks([0] + [i + 0.5 for i, state in enumerate(states)],['empty']+states)
    cbar.set_ylim(0,len(states))
    cbar.set_position([0.95, 0.2, 0.02,0.25])
    cbar.set_title('Class')
    plt.suptitle(suptitle,y=1)
    return {'fig':fig,'ax':ax,'cbar':cbar,'state_to_index':state_to_index, 'plate_indices':plate_indices}



def plot_plate(coordinate_values: dict,
               log: bool=True,
               vmin: float=None,
               vmax: float =None,
               cmap_name:str ='viridis',
               usenorm: bool=True, # Use normlizer (disable when using a custom colormap with discrete values
               cmap=None,
               colorbarargs={},
               returncb=False,
               class_colors: dict = None # When supplied no colormap is used, instead this mapping is used and a legend created
               ):



    coordinate_values  = {
        kwgs[:2]:value
        for kwgs, value in coordinate_values.items()
    }

    fig, ax = plt.subplots()
    n_rows = 16
    n_cols = 24
    if class_colors is None:
        if cmap is None:
            cmap = matplotlib.cm.get_cmap(cmap_name)

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
    ###########

    if class_colors is None:
        if vmax is None:
            vmax = np.percentile( list(coordinate_values.values()), 98)
            if log:
                vmax = np.power(10,np.ceil(np.log10(vmax)))

        if usenorm:
            if log:
                norm = matplotlib.colors.LogNorm(vmin=1 if vmin is None else vmin, vmax=vmax)
            else:
                norm = matplotlib.colors.Normalize(vmin=0 if vmin is None else vmin, vmax=vmax)

    for row,col in product(range(n_rows), range(n_cols)) :

        #if (y,x) in coordinate_values:
        #    print(np.clip(coordinate_values.get((y,x))/vmax,0,1))
        #print(None if (y,x) not in coordinate_values else cmap( np.clip(coordinate_values.get((y,x))/vmax,0,1)))

        if class_colors is not None:
            fc = class_colors[coordinate_values.get( (row,col), np.nan)]
        elif usenorm is False:
            fc = cmap(coordinate_values.get( (row,col), np.nan))
        else:
            fc = cmap(norm(coordinate_values.get( (row,col), np.nan)))

        ax.add_patch( Circle( (col,n_rows-row-1),
                             radius=0.45,
                             fill= (True if (row,col) not in coordinate_values else True),
                             fc= fc))

    ax.set_ylim(-1, n_rows)
    ax.set_xlim(-1, n_cols)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(np.arange(1,n_cols+1))

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([string.ascii_uppercase[n_rows-i-1] for i in range(n_rows)])
    ax.xaxis.tick_top()

    ax.grid()


    ax.tick_params(axis='y', which='both',length=3, pad=6)

    for label in ax.get_yticklabels():
        label.set_horizontalalignment('center')

    #norm = mpl.colors.Normalize(vmin=0, vmax=vmax)

    if class_colors is not None:
        patches = []
        for class_name, color in class_colors.items():
            patches.append( mpatches.Patch(color=color, label=class_name) )

        cb = legend = ax.legend(handles=patches,loc="upper left", bbox_to_anchor=(1,1))
        cax= None
    else:
        cax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,
                                        norm=norm if usenorm else None,
                                        orientation='vertical',**colorbarargs)

        cb.outline.set_visible(False)
    if returncb:
        return fig, ax, cax, cb
    else:
        return fig, ax, cax
