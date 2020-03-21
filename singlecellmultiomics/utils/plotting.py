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


    def get_relative_widths(self):
        return [self.lengths[contig]/self.total_bp for contig in self.contigs]

    def reset_axis(self, contig):
        ax = self[contig]
        ax.clear()

        #plt.axis('on')
        #ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_yticks([],[])

        #ax.set_aspect('equal')

        ax.set_xlabel(contig.replace('chr',''))
        ax.set_xlim(0,self.lengths[contig])

    def get_figure(self, figsize=(20,1)):

        widths = self.get_relative_widths()

        gs_kw = dict(width_ratios=widths)
        figure = plt.figure(figsize =figsize)
        self.gridspec = gridspec.GridSpec(1, len(widths), figure=figure, wspace=0.1, width_ratios=widths)
        #return plt.subplots(ncols=len(widths), nrows=1, constrained_layout=True,
        #        gridspec_kw=gs_kw, figsize=(20,4), sharey=True)
        self.axis = {}
        prev_ax = None
        for i,contig in enumerate(self.contigs):
           # i = i + 1 # grid spec indexes from 0

            ax = plt.subplot(self.gridspec[i], sharey=prev_ax)
            self.axis[contig] = ax

            self.reset_axis(contig)
            prev_ax=ax
        sns.despine()
        figure.canvas.draw()
        return figure

    def __getitem__(self, contig):
        return self.axis[contig]
