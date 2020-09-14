#!/usr/bin/env python3

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from pysamiterators import CachedFasta, MatePairIterator

# Molecule modules:
from singlecellmultiomics.molecule import TranscriptMolecule, MoleculeIterator
from singlecellmultiomics.fragment import SingleEndTranscriptFragment
from singlecellmultiomics.features import FeatureContainer

# Conversion modules:
from singlecellmultiomics.variants.substitutions import conversion_dict_stranded
from singlecellmultiomics.variants import substitution_plot, vcf_to_position_set
from singlecellmultiomics.utils import reverse_complement, complement
from collections import defaultdict, Counter
from singlecellmultiomics.utils import is_main_chromosome
from singlecellmultiomics.bamProcessing import sorted_bam_file, merge_bams


from scipy import stats
from multiprocessing import Pool
import os

import argparse

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def substitution_plot_stranded(pattern_counts: dict,
                      figsize: tuple = (12, 4),
                      conversion_colors: tuple = ('b', 'k', 'r', 'grey', 'g', 'pink','b','k','r','k','w','g'),
                      ylabel: str = '# conversions per molecule',
                      add_main_group_labels: bool = True,
                               ax=None,fig=None,
                      **plot_args
                      ):
    """
    Create 3bp substitution plot

    Args:
        pattern_counts(OrderedDict) : Dictionary containing the substitutions to plot.
            Use variants.vcf_to_variant_contexts to create it.
            Format:
            ```OrderedDict([(('ACA', 'A'), 0),
                 (('ACC', 'A'), 1),
                 (('ACG', 'A'), 0),
                 ...
                 (('TTG', 'G'), 0),
                 (('TTT', 'G'), 0)])```

        figsize(tuple) : size of the figure to create

        conversion_colors(tuple) : colors to use for the conversion groups

        ylabel(str) : y axis label

        add_main_group_labels(bool) : Add conversion group labels to top of plot

         **plot_args : Additional argument to pass to .plot()

    Returns
        fig : handle to the figure
        ax : handle to the axis

    Example:
        >>> from singlecellmultiomics.variants import vcf_to_variant_contexts, substitution_plot
        >>> import matplotlib.pyplot as plt
        >>> pobs = vcf_to_variant_contexts('variants.vcf.gz', 'reference.fasta')
        >>> for sample, conversions in pobs.items():
        >>>     fig, ax = substitution_plot(conversions)
        >>>     ax.set_title(sample)
        >>>     plt.show()

    """

    conversions_single_nuc = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')

    # Colors for the conversion groups:
    color_d = dict(zip(conversions_single_nuc, conversion_colors))
    colors = [color_d.get(f'{context[1]}{to}') for context, to in pattern_counts.keys()]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    substitution_dataframe = pd.DataFrame(pattern_counts.values(), index=list(pattern_counts.keys())).T
    substitution_dataframe.plot(kind='bar', color=colors, legend=False, width=1.0, ax=ax, edgecolor='k', **plot_args)
    offset = (1 / len(pattern_counts)) * 0.5  # Amount of distance for a half bar

    # Add 3bp context ticks:
    ax.set_xticks(np.linspace(-0.5 + offset, 0.5 - offset, len(pattern_counts)))
    ax.set_xticklabels( [context for context, to in pattern_counts.keys()], rotation=90, size=5)
    ax.set_ylabel(ylabel)
    ax.set_xlim((-0.5, 0.5))
    sns.despine()
    if add_main_group_labels:
        for i, (u, v) in enumerate(conversions_single_nuc):
            ax.text(  # position text relative to Axes
                (i + 0.5) / len(conversions_single_nuc), 1.0, f'{u}>{v}', fontsize=8,
                ha='center', va='top',
                transform=ax.transAxes,bbox=dict(facecolor='white', alpha=1,lw=0)
            )

    return fig, ax



if __name__=='__main__':

    argparser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      description='Assign molecules, set sample tags, set alleles')
    argparser.add_argument('bamin', type=str, help='Input BAM file')


    argparser.add_argument('-o', type=str, help="output bam file", required=True)


    argparser.add_argument('-reference', type=str, help="Reference_path", required=True)

    argparser.add_argument('-known', type=str, help="Known variants (vcf)", required=True)

    argparser.add_argument('-exons', type=str, help="Known variants (vcf)", required=True)

    argparser.add_argument('-introns', type=str, help="Known variants (vcf)", required=True)

    args = argparser.parse_args()


    single_cell_bam_path = args.bamin
    reference_path = args.reference

    # Load known variation, to ignore for mut-spectrum
    known_vcf_path = args.known

    # Paths to gene models
    exons_gtf_path = args.exons
    introns_gtf_path = args.introns

    # Write a tagged bam file to this path:
    tagged_output_path = args.o

    #####



    def obtain_conversions(contig : str):

        """ Create conversion dictionary for the suppled contig

        Args:
            contig (str)

        Returns:
            conversions_per_library (defaultdict( conversion_dict_stranded ) ) : Per library conversion dictionary
            n_molecules_per_library (Counter) : observed molecules per library
            contig(str) : the contig passed to the method
            temp_bam_path(str) : path to tagged bam file, tagged with gene annotations and 4su mutation count

        """

        conversions_per_library = defaultdict( conversion_dict_stranded )
        n_molecules_per_library = Counter()

        from singlecellmultiomics.molecule import might_be_variant

        # Create temp directory to write tagged bam file to:
        temp_dir = 'scmo_temp'
        temp_bam_path = f'{temp_dir}/{contig}.bam'
        if not os.path.exists('scmo_temp'):
            try:
                os.makedirs()
            except Exception as e:
                pass

        # Load gene annotations for the selected contig:
        transcriptome_features = FeatureContainer()
        transcriptome_features.loadGTF(
            path=exons_gtf_path,
            select_feature_type=['exon'],
            identifierFields=(
                'exon_id',
                'gene_id'),
            store_all=True,
            contig=contig,
            head=None)

        transcriptome_features.loadGTF(
            path=introns_gtf_path,
            select_feature_type=['intron'],
            identifierFields=['transcript_id'],
            store_all=True,
            contig=contig,
            head=None)


        try:
            with pysam.AlignmentFile(single_cell_bam_path) as alignments, \
                 pysam.VariantFile(known_vcf_path) as known, \
                 sorted_bam_file(temp_bam_path, origin_bam=single_cell_bam_path) as out, \
                 pysam.FastaFile(reference_path) as reference_handle:

                # Cache the sequence of the contig: (faster)
                reference = CachedFasta(reference_handle)

                for n_molecules, molecule in enumerate(MoleculeIterator(alignments,
                                                 TranscriptMolecule,
                                                 SingleEndTranscriptFragment,
                                                 fragment_class_args = {
                                                    'stranded':True,
                                                    'features':transcriptome_features
                                                },
                                                 molecule_class_args={
                                                     'reference':reference,
                                                     'features':transcriptome_features,
                                                     'auto_set_intron_exon_features': True
                                                 }, contig=contig
                                                 )):
                    # Read out mut spectrum
                    consensus = molecule.get_consensus()
                    n_molecules_per_library[molecule.library] += 1

                    n_4su_mutations = 0
                    n_4su_contexts = 0

                    for (chrom,pos), base in consensus.items():
                        context = reference.fetch(chrom, pos-1, pos+2).upper()


                        if ( (context[1]=='A' and  not molecule.strand) or (context[1]=='T' and molecule.strand) ) :
                            n_4su_contexts+=1


                        # Check if the base matches or the refence contains N's
                        if context[1]==base or 'N' in context or len(context)!=3:
                            continue

                        # Ignore germline variants:
                        if might_be_variant(chrom, pos,  known):
                            continue

                        if not molecule.strand: # reverse template
                            context = reverse_complement(context)
                            base = complement(base)

                        # Count 4SU specific mutations, and write to molecule later
                        if context[1]=='T' and base=='C':
                            n_4su_mutations+=1

                        conversions_per_library[molecule.library][(context, base)] += 1

                    # Write 4su modification to molecule
                    molecule.set_meta('4S',n_4su_mutations)
                    molecule.set_meta('4c',n_4su_contexts)
                    molecule.write_tags()
                    # Write tagged molecule to output file
                    molecule.write_pysam(out)

        except KeyboardInterrupt:
            # This allows you to cancel the analysis (CTRL+C) and get the current result
            pass

        return conversions_per_library, n_molecules_per_library, contig, temp_bam_path


####



    n_molecules_per_library = Counter()
    with Pool() as workers:
        conversions_per_library = defaultdict( conversion_dict_stranded ) # library : (context, query) : obs (int)

        # Obtain all contigs from the input bam file, exclude scaffolds:
        with pysam.AlignmentFile(single_cell_bam_path) as alignments:
            contigs = [contig for contig in alignments.references if is_main_chromosome(contig) and contig not in ['MT','Y'] ]

        # Run conversion detection on all contigs in parallel:
        tagged_bams = []
        for conversions_for_contig, \
            n_molecules_for_contig_per_lib, \
            contig, \
            temp_tagged_bam in workers.imap_unordered(obtain_conversions, contigs):

            # Merge conversion dictionary:
            for library, library_convs in conversions_for_contig.items():
                for context, observations in library_convs.items():
                    conversions_per_library[library][context] += observations

            n_molecules_per_library+=n_molecules_for_contig_per_lib

            print(f'finished {contig}  ', end='\r')

            tagged_bams.append(temp_tagged_bam)

        # Merge:
        print(f'Merging   ', end='\r')

        merge_bams(tagged_bams, tagged_output_path)

        # Normalize observed counts to the amount of molecules we saw:
        for library, library_convs in conversions_per_library.items():
            for context, observations in library_convs.items():
                library_convs[context] = observations / n_molecules_per_library[library]



    try:
        fig, axes = plt.subplots(len(conversions_per_library),1, figsize=(16,22), sharey=True )

        for ax, (library, conversions) in zip(axes,conversions_per_library.items()):

            substitution_plot_stranded(conversions,fig=fig, ax=ax,ylabel='conversions seen per molecule')

            ax.set_axisbelow(True)
            ax.grid(axis='y')
            ax.set_title(f'{library}, {n_molecules_per_library[library]} molecules')

        fig.tight_layout(pad=3.0)
        plt.savefig(f'conversions.png')
    except Exception as e:
        pass
