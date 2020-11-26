#!/usr/bin/env python

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
import pickle
import gzip
from uuid import uuid4

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
      description='Assign molecules')
    argparser.add_argument('bamin', type=str, help='Input BAM file')


    argparser.add_argument('-o', type=str, help="output bam file (.bam)", required=True)

    argparser.add_argument('-reference', type=str, help="Reference_path (.fasta)", required=True)

    argparser.add_argument('-known', type=str, help="Known variants (vcf)", required=True)

    argparser.add_argument('-exons', type=str, help="exons (gtf.gz)", required=True)

    argparser.add_argument('-introns', type=str, help="introns (gtf.gz)", required=True)

    argparser.add_argument('--R2_based', help="The input is only R2 sequences, the molcule mapping direction will be inverted", action='store_true')

    argparser.add_argument('-temp_dir', type=str, help="scmo_temp", default=str(uuid4()))

    argparser.add_argument('-tagthreads', type=int, help="Amount of threads used (int)", required=True)



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
        temp_dir = args.temp_dir
        temp_bam_path = f'{temp_dir}/{contig}.bam'
        if not os.path.exists(temp_dir):
            try:
                os.makedirs(temp_dir)
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



        colormap = plt.get_cmap('RdYlBu_r')
        colormap.set_bad((0,0,0))

        read_groups = {}
        try:
            with pysam.AlignmentFile(single_cell_bam_path, threads=4) as alignments, \
                 pysam.VariantFile(known_vcf_path) as known, \
                 sorted_bam_file(temp_bam_path, origin_bam=single_cell_bam_path, read_groups=read_groups, fast_compression=True) as out, \
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
                    if args.R2_based:
                        molecule.strand =  not molecule.strand # Invert becayse its R2 based.
                    n_molecules_per_library[molecule.library] += 1

                    n_4su_mutations = 0
                    n_4su_contexts = 0

                    for (chrom,pos), base in consensus.items():
                        context = reference.fetch(chrom, pos-1, pos+2).upper()
                        if len(context)!=3:
                            continue

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
                    # Set read color based on conversion rate:

                    try:
                        # The max color value will be 10% modification rate
                        cfloat = colormap( np.clip( 10*(n_4su_mutations/n_4su_contexts),0,1) )[:3]
                    except Exception as e:
                        cfloat = colormap._rgba_bad[:3]
                    molecule.set_meta('YC', '%s,%s,%s' % tuple((int(x * 255) for x in cfloat)))


                    molecule.set_meta('4c',n_4su_contexts)
                    molecule.write_tags()

                    for fragment in molecule:
                        rgid = fragment.get_read_group()
                        if not rgid in read_groups:
                            read_groups[rgid] = fragment.get_read_group(True)[1]

                    # Write tagged molecule to output file
                    molecule.write_pysam(out)

        except KeyboardInterrupt:
            # This allows you to cancel the analysis (CTRL+C) and get the current result
            pass

        return conversions_per_library, n_molecules_per_library, contig, temp_bam_path


    n_molecules_per_library = Counter()
    with Pool(args.tagthreads) as workers:
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
        fig, axes = plt.subplots(len(conversions_per_library),1, figsize=(16,4*(len(conversions_per_library))), sharey=True )
        if len(conversions_per_library)==1:
            axes = [axes]
        for ax, (library, conversions) in zip(axes,conversions_per_library.items()):
            # Export;:
            substitution_dataframe = pd.DataFrame(conversions.values(), index=list(conversions.keys())).T
            substitution_dataframe.to_csv(tagged_output_path.replace('.bam',f'{library}_conversions.csv'))
            substitution_plot_stranded(conversions,fig=fig, ax=ax,ylabel='conversions seen per molecule')

            ax.set_axisbelow(True)
            ax.grid(axis='y')
            ax.set_title(f'{library}, {n_molecules_per_library[library]} molecules')

        fig.tight_layout(pad=3.0)
        plt.savefig(tagged_output_path.replace('.bam','conversions.png'))
    except Exception as e:
        print(e)




    # Count amount of 4sU conversions per cell, per gene
    def listionary():
        return defaultdict(list)

    expression_per_cell_per_gene = defaultdict(Counter) # gene -> cell -> obs
    four_su_per_cell_per_gene = defaultdict(listionary ) # cell -> gene -> [] 4_su observation counts per molecule
    four_su_per_gene_per_cell = defaultdict(listionary ) # gene -> cell -> [] 4_su observation counts per molecule

    with pysam.AlignmentFile(tagged_output_path) as reads:
        for R1,R2 in MatePairIterator(reads):

            for read in (R1,R2): # Count every fragment only once by selecting one of the two reads.
                if read is not None:
                    break

            if read.has_tag('gn'):
                gene = read.get_tag('gn')
            elif read.has_tag('GN'):
                gene = read.get_tag('GN')
            else:
                continue

            if read.is_duplicate:
                continue

            cell = read.get_tag('SM')
            foursu = read.get_tag('4S')
            foursu_contexts = read.get_tag('4c')
            library = read.get_tag('LY')
            cell = cell.split('_')[1] # Remove library part
            expression_per_cell_per_gene[gene][(library,cell)] += 1
            if foursu_contexts>0:
                four_su_per_gene_per_cell[gene][(library,cell)].append(foursu/foursu_contexts)
                four_su_per_cell_per_gene[(library,cell)][gene].append(foursu/foursu_contexts)
            assert not (foursu>0 and foursu_contexts==0)

    # Store these dictionaries to disk
    with gzip.open( tagged_output_path.replace('.bam','4sU_per_gene_per_cell.dict.pickle.gz'),'wb' ) as o:
        pickle.dump(four_su_per_gene_per_cell, o)
    with gzip.open( tagged_output_path.replace('.bam','4sU_per_cell_per_gene.dict.pickle.gz'),'wb' ) as o:
        pickle.dump(four_su_per_cell_per_gene, o)

    with gzip.open( tagged_output_path.replace('.bam','expression_per_cell_per_gene.pickle.gz'),'wb' ) as o:
        pickle.dump(expression_per_cell_per_gene, o)

    four_su_per_gene_per_cell_mean = defaultdict(dict)
    four_su_per_gene_per_cell_total= defaultdict(dict)
    for gene  in four_su_per_gene_per_cell:
        for cell, fsu_obs in four_su_per_gene_per_cell[gene].items():
            four_su_per_gene_per_cell_mean[gene][cell] = np.mean(fsu_obs)
            four_su_per_gene_per_cell_total[gene][cell] = np.sum( np.array(fsu_obs)>0 )
    four_su_per_gene_per_cell_mean = pd.DataFrame(four_su_per_gene_per_cell_mean).T
    four_su_per_gene_per_cell_total = pd.DataFrame(four_su_per_gene_per_cell_total).T

    four_su_per_gene_per_cell_mean.to_csv(tagged_output_path.replace('.bam','4sU_labeled_ratio.csv.gz'))

    expression_matrix = pd.DataFrame(four_su_per_gene_per_cell).T.fillna(0)
    libraries = expression_matrix.columns.get_level_values(0).unique()

    ############
    fig, ax = plt.subplots(figsize=(7,7))

    min_molecules = 100

    conversion_ratios = {} # cell->gene->ratio

    for library in sorted(list(libraries)):

        if not '4s' in library:
            continue

        cell_efficiencies = {}
        cell_molecules = Counter()



        for cell, genes in  four_su_per_cell_per_gene.items():
            target_cell_name = cell
            if cell[0]!=library:
                continue
            if '100cells' in library:
                target_cell_name = 'bulk'

            conversions_total = []
            for gene, conversions in genes.items():
                conversions_total+= conversions
                cell_molecules[target_cell_name]+=len(conversions)
            cell_efficiencies[target_cell_name] = np.mean(conversions_total)*100



        selected_cells = [cell for cell in cell_efficiencies if cell_molecules[cell]>min_molecules]

        cell_efficiencies = {cell:cell_efficiencies[cell] for cell in selected_cells}

        scatter = plt.scatter( [cell_molecules[cell] for cell in selected_cells],
                              cell_efficiencies.values(),
                              label=library, s=2 )
        plt.scatter(
              np.median([cell_molecules[cell] for cell in selected_cells]),
               np.median( list(cell_efficiencies.values())),
            facecolors = 'k',
            s=250,
            marker='+',edgecolors='black', lw=3
        )
        plt.scatter(
              np.median([cell_molecules[cell] for cell in cell_efficiencies]),
               np.median( list(cell_efficiencies.values())),
            facecolors = scatter.get_facecolors(),
            s=250,
            marker='+',edgecolors='black', lw=1
        )



    plt.ylabel('4sU conversion rate (%)')
    plt.xlabel('total molecules')
    plt.xscale('log')
    ax.set_axisbelow(True)
    ax.grid()
    plt.legend( bbox_to_anchor=(0.6, 1))
    sns.despine()
    plt.title('4sU conversion rate per cell')
    plt.savefig(tagged_output_path.replace('.bam','conversion_rate.png'), dpi=200)


    ##########



    fig, axes = plt.subplots(6,4,figsize=(13,17), squeeze=True)
    axes = axes.flatten()
    axes_index = 0



    for gene in expression_matrix.mean(1).sort_values()[-100:].index:

        fraction_hits = defaultdict(list)
        labeled = defaultdict(list)
        total = defaultdict(list)


        if gene=='MALAT1':
            continue

        for (library, cell), labeled_4su_fraction in four_su_per_gene_per_cell[gene].items():

            #if not '4sU' in library and not 'LIVE' in library:
            #    continue

            if '4sU' in library:
                library = '4sU'
            else:
                library= 'unlabeled'

            fraction_hits[library] += labeled_4su_fraction
            labeled[library].append( sum([ l>0 for l in labeled_4su_fraction]) )
            total[library].append(len(labeled_4su_fraction))


        try:
            max_x = max( ( max(total[library]) for library in total))


            slope, intercept, r_value, p_value, std_err = stats.linregress(total['4sU'],labeled['4sU'])
            if slope<0.001 or p_value>0.05 or np.isnan(p_value):
                continue

            slope, intercept, r_value, p_value, std_err = stats.linregress(total['unlabeled'],labeled['unlabeled'])
            if p_value>0.05 or np.isnan(p_value) :
                continue
        except Exception as e:
            continue
        ax = axes[axes_index]
        axes_index+=1


        for library in total:
            slope, intercept, r_value, p_value, std_err = stats.linregress(total[library],labeled[library])

            #max_x = max(total[library])
            ax.plot([0,max_x],[intercept,max_x*slope + intercept],c='red' if '4sU' in library else 'k' )




        for library in total:
            ax.scatter(total[library],labeled[library], label=library , s=10, alpha=0.5, c='red' if '4sU' in library else 'k' )



        slope, intercept, r_value, p_value, std_err = stats.linregress(total['4sU'],labeled['4sU'])
        ax.legend()
        ax.set_xlabel('total molecules')
        ax.set_ylabel('4sU labeled molecules')
        ax.set_title(f'{gene}\nslope:{slope:.2f}')
        sns.despine()

        if axes_index>=len(axes):
            break

        fig.tight_layout(pad=1.0)
        plt.savefig( (tagged_output_path.replace('.bam','slopes.png')) , dpi=200)
