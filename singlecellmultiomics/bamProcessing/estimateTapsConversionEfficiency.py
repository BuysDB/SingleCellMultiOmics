#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from singlecellmultiomics.molecule import MoleculeIterator, CHICNLAMolecule, TAPSNlaIIIMolecule,TAPSCHICMolecule,TAPS
from singlecellmultiomics.fragment import CHICFragment, NlaIIIFragment
import pysam
from pysamiterators import CachedFasta
from singlecellmultiomics.variants.substitutions import conversion_dict_stranded
from collections import defaultdict
from singlecellmultiomics.utils import reverse_complement, complement
from glob import glob
from multiprocessing import Pool
from singlecellmultiomics.bamProcessing.bamFunctions import get_reference_path_from_bam
from collections import Counter
import pandas as pd
import matplotlib as mpl
from singlecellmultiomics.utils.sequtils import phred_to_prob, prob_to_phred
import seaborn as sns
import argparse

def update_mutation_dict(molecule,reference, conversions_per_library, context_obs):

    consensus = molecule.get_consensus(dove_safe=True,
                                       min_phred_score=22,
                                       skip_first_n_cycles_R1=10,
                                       skip_last_n_cycles_R1=20,
                                       skip_first_n_cycles_R2=10,
                                       skip_last_n_cycles_R2=20,
                                       dove_R2_distance=15,
                                       dove_R1_distance=15


                                      )

    nm = 0


    contexts_to_add = []


    for (chrom,pos), base in consensus.items():
        context = reference.fetch(chrom, pos-1, pos+2).upper()

        if len(context)!=3:
            continue

        # Check if the base matches or the refence contains N's
        if 'N' in context or len(context)!=3:
            continue

        # Ignore germline variants:
        #if might_be_variant(chrom, pos,  known):
        #    continue

        if not molecule.strand: # reverse template
            context = reverse_complement(context)
            base = complement(base)

        if context[1]!='C' and context[1]!=base:
            nm+=1

        contexts_to_add.append((context,base))


    if nm>5:
        nm=5

    k = tuple((*molecule.sample.rsplit('_',2), nm))
    for (context, base) in contexts_to_add:

        context_obs[ k ][context] += 1
        try:
            conversions_per_library[k][(context, base)] += 1
        except:
            pass


def get_conversion_counts(args):


    taps = TAPS()

    conversions_per_library = defaultdict( conversion_dict_stranded )
    context_obs = defaultdict( Counter )

    bam,refpath,method,every_fragment_as_molecule, spikein_name  = args


    if method=='nla':
        fragment_class=NlaIIIFragment
        molecule_class=TAPSNlaIIIMolecule
    else:
        fragment_class=CHICFragment
        molecule_class=TAPSCHICMolecule

    with pysam.FastaFile(refpath) as ref:
        reference = CachedFasta(ref)


        print(f'Processing {bam}')

        with pysam.AlignmentFile(bam, threads=8) as al:

            for molecule in MoleculeIterator(
                al,
                fragment_class=fragment_class,
                molecule_class=molecule_class,
                molecule_class_args={
                     'reference':reference,
                     'taps':taps,
                    'taps_strand':'R'
                },
                every_fragment_as_molecule=every_fragment_as_molecule,
                fragment_class_args={},
                contig = spikein_name
            ):
                update_mutation_dict(molecule, reference ,conversions_per_library, context_obs)

    return conversions_per_library, context_obs


def generate_taps_conversion_stats(bams, reference_path, prefix, method, every_fragment_as_molecule, spikein_name, n_threads=None):
    if reference_path is None:
        reference_path = get_reference_path_from_bam(bams[0])

    print(f'Reference at {reference_path}')
    if reference_path is None:
        raise ValueError('Please supply a reference fasta file')

    conversions_per_library = defaultdict( conversion_dict_stranded )
    context_obs = defaultdict( Counter )

    with Pool(n_threads) as workers:

        for cl, co in workers.imap(get_conversion_counts, [(bam, reference_path, method, every_fragment_as_molecule, spikein_name) for bam in bams] ):

            for lib, obs in cl.items():
                for k,v in obs.items():
                    conversions_per_library[lib][k] +=v

            for lib, obs in co.items():
                for k,v in obs.items():
                    context_obs[lib][k] += v


    qf = pd.DataFrame(context_obs)
    qf.to_csv(f'{prefix}_conversions_counts_raw_lambda.csv')
    ###
    indices = []
    for lib, qqf in qf.groupby(level=0,axis=1):

        ser = qqf.sum(level=(0,1),axis=1).sum().sort_values(ascending=False)
        ser = ser[ser>5000][::-1]
        indices += list(ser.index)


    ###

    normed_conversions_per_library = defaultdict( conversion_dict_stranded )

    for INDEX in context_obs:
        for (context, base),obs in conversions_per_library[INDEX].items():
            try:
                normed_conversions_per_library[INDEX][(context,base)] = obs/ context_obs[INDEX][context]
            except Exception:
                pass

    df = pd.DataFrame(normed_conversions_per_library)
    df = df[ [INDEX for INDEX in df if (INDEX[0],  INDEX[1]) in indices] ]

    df = df.loc[ [(context, base)for context, base in df.index if context[1]=='C' and base=='T' and context.endswith('CG')]  ]
    df = df.T

    df.to_csv(f'{prefix}_conversions_lambda.csv')

    mpl.rcParams['figure.dpi'] = 300

    samples = []

    for (lib, cell, nm), row in df.iterrows():

        if nm!=0:
            continue

        for context, base in [('ACG', 'T'),
            ('CCG', 'T'),
            ('GCG', 'T'),
            ('TCG', 'T')]:

            r = {
            'lib':lib,
            'cell':cell,
            'nm':nm,
                #'plate': int(lib.split('-')[-1].replace('pl','')),
            'group':  f'{nm},{context},{cell}',
            'context': f'{context}>{base}',
            'conversion rate':row[context,base]
            }

            samples.append(r)

    plot_table = pd.DataFrame(samples)
    print(plot_table)

    ph = 22
    fig, ax = plt.subplots(figsize=(12,5))
    sns.boxplot(data=plot_table.sort_values('lib'),x='context', y='conversion rate',hue='lib',whis=6, ax=ax)

    #ax = sns.swarmplot(data=plot_table,x='nm', y='conversion rate',hue='lib',)

    plt.legend()
    plt.ylabel('Lambda Conversion rate')
    sns.despine()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5 ))

    plt.suptitle('Estimated TAPS conversion rate', y=1.05, fontsize=12)
    plt.title(f'Lambda spike-in, >{(1.0-(phred_to_prob(22)))*100 : .2f}% accuracy base calls', fontsize=10)
    plt.tight_layout()
    plt.savefig(f'{prefix}_conversion_rate_phred_{ph}.png', bbox_inches='tight')
    plt.close()


if __name__=='__main__':
    argparser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      description='Estimate the conversion efficiency of a TAPS converted file. ')
    argparser.add_argument('bams', type=str, nargs='+',help='Input bam files')
    argparser.add_argument('-o', type=str, help="output alias (Will be the prefix of the output files)", required=True)
    argparser.add_argument('-method', type=str, default='nla', help='Molecule class (nla or chic). Use chic when you are not sure or when another other protocol is used.')
    argparser.add_argument('--dedup', action='store_true',help='perform UMI deduplication and consensus calling. Do not use when the UMI\'s are (near) saturated')
    argparser.add_argument('-t', type=int, help='Amount of threads')
    argparser.add_argument('-spikein_name', type=str, help='Name of spikein contig',default='J02459.1')


    argparser.add_argument(
        '-ref',
        type=str,
        default=None,
        help="Path to reference fast (autodected if not supplied)")
    args = argparser.parse_args()


    generate_taps_conversion_stats(args.bams,
                                   args.ref,
                                   prefix=args.o,
                                   method=args.method,
                                   every_fragment_as_molecule=not args.dedup,
                                   spikein_name=args.spikein_name,
                                   n_threads=args.t)
