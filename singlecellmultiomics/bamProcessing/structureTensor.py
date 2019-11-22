#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import itertools
import random
import imageio
import numpy as np
import sys
import os
import argparse
import singlecellmultiomics.bamProcessing
import collections
import pysamiterators
import pandas as pd
import pysam
import singlecellmultiomics.fragment
import singlecellmultiomics.molecule
from singlecellmultiomics.molecule import MoleculeIterator
import singlecellmultiomics
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')


def empty_tensor_dict():
    d = {
        'insert': 0,
        'deleted': 0,
        'match': 0,
        'mismatch': 0,
        'clip': 0,
        'A': 0,
        'T': 0,
        'C': 0,
        'G': 0}
    for baseA, baseB in itertools.product('ACTG', repeat=2):
        d[f'{baseA}>{baseB}'] = 0
    return d


class TensorStackMolecule(singlecellmultiomics.molecule.Molecule):

    def __init__(self,
                 fragments=None,
                 cache_size=10_000,
                 reference=None,
                 min_max_mapping_quality=None,
                 # When all fragments have a mappin quality below this value
                 # the is_valid method will return False
                 allele_resolver=None,
                 **kwargs):
        singlecellmultiomics.molecule.Molecule.__init__(self, fragments,
                                                        cache_size,
                                                        reference,
                                                        min_max_mapping_quality,
                                                        # When all fragments
                                                        # have a mappin quality
                                                        # below this value the
                                                        # is_valid method will
                                                        # return False
                                                        allele_resolver,
                                                        **kwargs)

    def get_alignment_tensor(self,
                             max_reads,
                             window_radius=20,
                             centroid=None,
                             mask_centroid=False,
                             refence_backed=False,
                             skip_missing_reads=False
                             ):
        """ Obtain a tensor representation of the molecule alignment around the given centroid
        Args:
            max_reads (int) : maximum amount of reads returned in the tensor, this will be the amount of rows/4 of the returned feature matrix

            window_radius (int) : radius of bp around centroid

            centroid(int) : center of extracted window, when not specified the cut location of the molecule is used

            mask_centroid(bool) : when True, mask reference base at centroid with N

            refence_backed(bool) : when True the molecules reference is used to emit reference bases instead of the MD tag
        Returns:
            tensor_repr(np.array) : (4*window_radius*2*max_reads) dimensional feature matrix
        """
        reference = None
        if refence_backed:
            reference = self.reference
            if self.reference is None:
                raise ValueError(
                    "refence_backed set to True, but the molecule has no reference assigned. Assing one using pysam.FastaFile()")

        height = max_reads
        chromosome = self.chromosome
        if centroid is None:
            _, centroid, strand = self.get_cut_site()
        span_start = centroid - window_radius
        span_end = centroid + window_radius
        span_len = abs(span_start - span_end)
        base_content_table = np.zeros((height, span_len))
        base_mismatches_table = np.zeros((height, span_len))
        base_indel_table = np.zeros((height, span_len))
        base_qual_table = np.zeros((height, span_len))
        base_clip_table = np.zeros((height, span_len))
        pointer = 0

        mask = None
        if mask_centroid:
            mask = set((chromosome, centroid))

        fragments = []
        for frag in self.fragments:
            _, frag_start, frag_end = frag.span
            if frag_start > span_end or frag_end < span_start:
                continue

            fragments.append(frag)

        print(len(fragments))
        for frag in random.sample(fragments, max_reads):

            pointer = frag.write_tensor(
                chromosome,
                span_start,
                span_end,
                index_start=pointer,
                base_content_table=base_content_table,
                base_mismatches_table=base_mismatches_table,
                base_indel_table=base_indel_table,
                base_qual_table=base_qual_table,
                base_clip_table=base_clip_table,
                height=height,
                mask_reference_bases=mask,
                reference=reference,
                skip_missing_reads=skip_missing_reads)

        x = np.vstack(
            [
                base_content_table,
                base_mismatches_table,
                base_indel_table,
                base_qual_table,
                base_clip_table
            ])

        return x


class SingleEndTranscriptTensorable(singlecellmultiomics.fragment.Fragment):
    def __init__(self, reads, **kwargs):
        singlecellmultiomics.fragment.Fragment.__init__(self, reads, **kwargs)

    def get_site_location(self):
        return self.get_span()[1]

    def umi_eq(self, umi):
        return True

    def __eq__(self, other):
        if min(abs(self.span[1] -
                   other.span[1]), abs(self.span[2] -
                                       other.span[2])) > self.assignment_radius:
            return False
        return True


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Assign molecules, set sample tags, set alleles')

    argparser.add_argument('bamin', type=str)
    argparser.add_argument('-contig', type=str)
    argparser.add_argument('-ref', type=str)
    argparser.add_argument('-mem', type=int, default=10)
    argparser.add_argument('-time', type=int, default=10)
    argparser.add_argument('-radius', type=int, default=5)
    #argparser.add_argument('-maxreads',  type=int, default=100)
    argparser.add_argument('-min_coverage', type=int, default=50)
    argparser.add_argument('-mindiff', type=float, default=0.04)
    argparser.add_argument(
        '-sample_non_variable_per_contig',
        type=int,
        default=2,
        help='Select this amount of non-variable bases per contig')
    argparser.add_argument('-head', type=int)
    argparser.add_argument('-o', type=str, default='./tensors')
    args = argparser.parse_args()
    if not os.path.exists(args.o):
        try:
            os.makedirs(args.o)
        except Exception as e:
            pass

    alignments_positive = pysam.AlignmentFile(args.bamin)

    if args.contig is None:
        try:
            contigs_per_job = 50
            contigs_in_job = []
            for i, stat in enumerate(
                    alignments_positive.get_index_statistics()):
                if stat.mapped < args.min_coverage:
                    continue
                contig = stat.contig
                contigs_in_job.append(contig)
                if len(contigs_in_job) >= contigs_per_job:
                    arguments = " ".join(
                        [x for x in sys.argv]) + f" -contig {','.join(contigs_in_job)}"
                    job = f'STRUCT_{i}'
                    os.system(
                        f'submission.py --silent' +
                        f' -y --py36 -time {args.time} -t 1 -m {args.mem} -N {job} " {arguments};"')
                    contigs_in_job = []
            # eject remains
            if len(contigs_in_job) > 0:
                arguments = " ".join([x for x in sys.argv]) + \
                    f" -contig {','.join(contigs_in_job)}"
                job = f'STRUCT_{i}'
                os.system(
                    f'submission.py --silent' +
                    f' -y --py36 -time {args.time} -t 1 -m {args.mem} -N {job} " {arguments};"')
                contigs_in_job = []

        except KeyboardInterrupt:
            exit()
    else:
        for contig in args.contig.split(','):

            print(f"Processing {contig}")

            non_variant_enough_cov = []
            variable_locations = []  # (pos)
            obs = {}
            for column in alignments_positive.pileup(contig=contig):
                edits = collections.Counter()
                for pileupread in column.pileups:
                    if pileupread.is_del:
                        edits['del'] += 1
                    elif pileupread.is_refskip:
                        edits['insert'] += 1
                    else:
                        edits[pileupread.alignment.query_sequence[pileupread.query_position]] += 1

                if column.nsegments >= args.min_coverage:

                    if edits.most_common(1)[0][1] < (
                            1 - args.mindiff) * sum(edits.values()):
                        variable_locations.append(column.pos)
                        obs[(contig, column.pos)] = edits
                        print(f'{column.pos} is variable')
                    else:
                        non_variant_enough_cov.append(column.pos)

            print('Calculating conversion tensor')
            conversions = collections.defaultdict(empty_tensor_dict)

            reference_bases = {}
            for ii, read in enumerate(alignments_positive.fetch(contig)):

                last_ref_pos = None
                clipped = 0
                for query_pos, ref_pos, ref_base in read.get_aligned_pairs(
                        with_seq=True):
                    if ref_pos is None:  # clip or skip
                        if last_ref_pos is not None:  # otherwise we don't know our starting coordinate
                            conversions[(contig, last_ref_pos)]['insert'] += 1
                        else:
                            clipped += 1
                    elif query_pos is None:
                        conversions[(contig, ref_pos)]['deleted'] += 1
                    else:
                        query_base = read.seq[query_pos]

                        if query_base != 'N':
                            conversions[(contig, ref_pos)][query_base] += 1
                            if ref_base != 'N':
                                conversions[(
                                    contig, ref_pos)][f'{ref_base.upper()}>{query_base}'] += 1
                        if query_base == ref_base.upper():
                            conversions[(contig, ref_pos)]['match'] += 1
                        else:
                            conversions[(contig, ref_pos)]['mismatch'] += 1

                    if ref_base is not None and ref_pos is not None:
                        reference_bases[(contig, ref_pos)] = ref_base
                    if ref_pos is not None and clipped > 0:
                        conversions[(contig, ref_pos)]['clip'] = 1
                        clipped = 0

                    last_ref_pos = ref_pos

                if args.head is not None and ii >= args.head - 1:
                    print('Head was supplied, stopped loading more reads')
                    break

            df = pd.DataFrame(conversions).fillna(0).T

            df['het2'] = [df.loc[index, ['A', 'T', 'C', 'G']
                                 ].T.nsmallest(2).sum() for index in df.index]
            df['het3'] = [df.loc[index, ['A', 'T', 'C', 'G']
                                 ].T.nsmallest(3).sum() for index in df.index]

            pd.DataFrame({'ref': reference_bases}).to_csv(
                f'{args.o}/{contig}.reference_bases.pickle.gz')
            df.to_pickle(f'{args.o}/{contig}.pickle.gz')
            df.to_csv(f'{args.o}/{contig}.feature_matrix.csv')

            X = []
            positions = []
            features_to_plot = [
                'insert',
                'deleted',
                'match',
                'mismatch',
                'clip']
            for location in variable_locations:
                radius = args.radius
                features = df.loc[contig].loc[location -
                                              radius:location + radius]
                features['is_variable'] = 1
                features[features_to_plot].plot()
                plt.savefig(f'{args.o}/{contig}.{location}.png')
                X.append(features.values.flatten())
                positions.append((contig, location))
                features.to_csv(f'{args.o}/{contig}.{location}.features.csv')

            # Select non-variable sites
            if len(non_variant_enough_cov) > 0:
                for location in random.sample(non_variant_enough_cov, min(
                        args.sample_non_variable_per_contig, len(non_variant_enough_cov))):
                    radius = args.radius
                    features = df.loc[contig].loc[location -
                                                  radius:location + radius]
                    features['is_variable'] = 0
                    features[features_to_plot].plot()
                    plt.savefig(f'{args.o}/{contig}.{location}.png')
                    X.append(features.values.flatten())
                    positions.append((contig, location))
                    features.to_csv(
                        f'{args.o}/{contig}.{location}.features.csv')

            pd.DataFrame(positions).to_csv(f'{args.o}/{contig}.positions.csv')
            pd.DataFrame(X).to_csv(f'{args.o}/{contig}.X.csv')

            """
            # Pool all reads in a single contig into one "molecule"
            print("Loading reads into memory...")
            molecule = TensorStackMolecule(reference=reference)

            for ii,R2 in enumerate(alignments_positive.fetch(contig=args.contig)):
                if R2.is_duplicate or R2.is_supplementary:
                    continue
                fragment = SingleEndTranscriptTensorable( (None, R2) )
                molecule._add_fragment(fragment)

                if ii%100_000==0:
                    print(f'Loaded {ii} reads', end='\r')
                if args.head is not None and ii>=args.head-1:
                    print('Head was supplied, stopped loading more reads')
                    break
            print('Finished loading reads to memory')
            pd.DataFrame(obs).to_csv(f'{args.o}/{contig}-SELECTED.csv' )
            for position in range(molecule.spanStart, molecule.spanEnd):
                if not position in variable_locations:
                    continue
                print(f'Writing tensor for {position}')
                x = molecule.get_alignment_tensor(centroid=position,window_radius=args.radius,max_reads=args.maxreads,skip_missing_reads=True, refence_backed=True)

                # write the tensor to a file
                pd.DataFrame(x).to_csv(f'{args.o}/{contig}-{position}.csv' )
                pd.DataFrame(x).to_pickle(f'{args.o}/{contig}-{position}.pickle.gz' )
                imageio.imwrite(f'{args.o}/{contig}-{position}.png', x)"""
