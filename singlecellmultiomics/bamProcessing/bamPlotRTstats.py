#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import singlecellmultiomics.features
import matplotlib.pyplot as plt
import numpy as np
import itertools
import pandas as pd
import importlib
import singlecellmultiomics.universalBamTagger.universalBamTagger as ut
import pysamiterators.iterators as pyts
import matplotlib.lines as mlines
import os
import sys
import pysam
import collections
import argparse
import gzip
import pickle
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def nlaIII_molecule_acceptance_function(molecule):
    first_read = molecule[0][0]
    if first_read.mapping_quality < 60:
        return False
    reject = False
    if first_read.has_tag('XA'):
        for alt_align in first_read.get_tag('XA').split(';'):
            if len(alt_align) == 0:  # Sometimes this tag is empty for some reason
                continue

            hchrom, hpos, hcigar, hflag = alt_align.split(',')
            if not hchrom.endswith('_alt'):
                reject = True
                break
    return reject


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Visualise feature density of a bam file. (Coverage around stop codons, start codons, genes etc)')
    argparser.add_argument('bamFile', type=str)
    argparser.add_argument(
        '-head',
        type=int,
        help='Process this many molecules')
    argparser.add_argument('-binSize', type=int, default=30)
    argparser.add_argument(
        '-maxfs',
        type=int,
        default=900,
        help='X axis limit of fragment size plot')
    argparser.add_argument('-maxOverseq', type=int, default=4)
    argparser.add_argument('--notstrict', action='store_true')
    argparser.add_argument('-o', type=str, default='RT_dist')

    args = argparser.parse_args()

    def dd():
        return collections.defaultdict(collections.Counter)

    fragment_distribution = collections.Counter()
    fragment_distribution_raw = collections.defaultdict(
        collections.Counter)  # overseq -> fragmentisze -> counts
    # Read size fragments
    fragment_distribution_raw_rf = collections.defaultdict(
        dd)  # lib -> overseq -> fragmentisze (span) -> counts
    gc_distribution = collections.Counter()
    gc_frag_distribution = collections.defaultdict(
        collections.Counter)  # fragment size -> observed gc/at+gc ratio
    # fragmentSize -> umi obs

    rt_frag_distribution = collections.defaultdict(collections.Counter)
    # fragment size -> amount of RT reactions -> count

    observed_cuts = collections.defaultdict()
    used = 0
    used_reads = 0
    gc_capture = False
    with pysam.AlignmentFile(args.bamFile) as a:
        for i, molecule in enumerate(
            ut.MoleculeIterator_OLD(
                a, umi_hamming_distance=1)):

            if not args.notstrict and not nlaIII_molecule_acceptance_function(
                    molecule):
                continue

            if args.head is not None and used > args.head:
                print('Stoppping, saw enough molecules')
                break

            rt_reactions = ut.molecule_to_random_primer_dict(molecule)
            amount_of_rt_reactions = len(rt_reactions)

            # this obtains the maximum fragment size:
            frag_chrom, frag_start, frag_end = pyts.getListSpanningCoordinates(
                [v for v in itertools.chain.from_iterable(molecule) if v is not None])

            # Obtain the fragment sizes of all RT reactions:
            rt_sizes = []
            for (rt_end, hexamer), fragment in rt_reactions.items():
                rt_chrom, rt_start, rt_end = pyts.getListSpanningCoordinates(
                    itertools.chain.from_iterable(fragment))
                rt_sizes.append([rt_end - rt_start])

            first_read = molecule[0][0]
            site = first_read.get_tag('DS')
            strand = first_read.get_tag('RS')
            library = first_read.get_tag('LY')
            # if gc_capture:
            #    sequence = reference.fetch(frag_chrom, frag_start, frag_end)
            #    gc = sequence.count('C')+ sequence.count('G')
        #        length = len(sequence)
            used += 1

            # if gc_capture:
            #    gc_distribution[gc/length] += 1
            #    gc_frag_distribution[fragment_size][gc/length] += 1
            #    fragment_distribution_raw[len(molecule)][fragment_size]+=1

            if len(rt_sizes) == 0:
                mean_rt_size = 0
            else:
                mean_rt_size = int(np.mean(rt_sizes))
            fragment_distribution_raw_rf[library][len(
                molecule)][mean_rt_size] += 1
            rt_frag_distribution[mean_rt_size][len(rt_reactions)] += 1
            used_reads += len(molecule)

    bin_size = args.binSize
    m_overseq = args.maxOverseq

    with gzip.open(f'{args.o}_raw_data.pickle.gz', 'wb') as fo:

        pickle.dump(fragment_distribution_raw_rf, fo)

    for library in fragment_distribution_raw_rf:
        fig, axes = plt.subplots(1, 1, figsize=(10, 7))
        ax = axes
        ax.set_title(
            f'Read fragment size distribution\n{used} molecules / {used_reads} fragments  analysed\n{library}')
        table = {}

        for overseq in range(1, m_overseq):
            try:
                # Rebin in 10bp bins:
                rebinned = collections.Counter()
                for f_size, obs in fragment_distribution_raw_rf[library][overseq].most_common(
                ):
                    rebinned[int(f_size / bin_size) * bin_size] += obs
                obs_dist_fsize = np.array(list(rebinned.keys()))
                obs_dist_freq = np.array(list(rebinned.values()))
                sorting_order = np.argsort(obs_dist_fsize)
                obs_dist_fsize = obs_dist_fsize[sorting_order]
                obs_dist_freq = obs_dist_freq[sorting_order]
                obs_dist_density = obs_dist_freq / np.sum(obs_dist_freq)

                for i, x in enumerate(obs_dist_fsize):
                    table[(overseq, x)] = {'obs_dist_freq': obs_dist_freq[i],
                                           'obs_dist_density': obs_dist_density[i]}

                ax.plot(
                    obs_dist_fsize, obs_dist_density, c=(
                        overseq / m_overseq, 0, 0), label=f'{overseq} read fragments / umi')
                ax.set_xlim(-10, args.maxfs)
                ax.legend()
                ax.set_xlabel('fragment size')
                ax.set_ylabel('density')
            except Exception as e:
                print(e)
                pass
        plt.tight_layout()
        plt.savefig(f'{args.o}_{library}.png')
        pd.DataFrame(table).to_csv(f'{args.o}_{library}.csv')
        try:
            plt.close()
        except Exception as e:
            pass

    # Export the table:
