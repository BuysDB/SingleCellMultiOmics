#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import gzip
import pickle
import matplotlib

import numpy as np

import pandas as pd
import singlecellmultiomics.features


def distance_to_hit(
        lookup_coordinate,
        hit_start,
        hit_end,
        hit_strand,
        lookup_strand=None,
        distance_from_begin=True):
    distance = None
    if distance_from_begin:
        if hit_strand == '+':
            distance = lookup_coordinate - hit_start
        else:
            distance = hit_start - lookup_coordinate
    else:
        if hit_strand == '+':
            distance = lookup_coordinate - hit_end
        else:
            distance = hit_end - lookup_coordinate

    return distance

# check if lookup coordinate falls within a feature:


def distance_to_feature_start(
        chromosome,
        lookup_coordinate,
        feature_container,
        lookup_strand=None):
    distance = None
    min_distance = None
    for hit in feature_container.findFeaturesAt(
            chromosome=chromosome,
            lookupCoordinate=lookup_coordinate,
            strand=lookup_strand):
        hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
        # calculate distance relative to selected point
        distance = distance_to_hit(
            lookup_coordinate,
            hit_start,
            hit_end,
            hit_strand,
            lookup_strand)
        if min_distance is None or abs(distance) < abs(min_distance):
            min_distance = distance
    if distance is None:
        for hit in feature_container.findNearestFeature(
                chromosome=chromosome,
                lookupCoordinate=lookup_coordinate,
                strand=lookup_strand):
            hit_start, hit_end, hit_id, hit_strand, hit_ids = hit

            distance = distance_to_hit(
                lookup_coordinate,
                hit_start,
                hit_end,
                hit_strand,
                lookup_strand)

            if min_distance is None or abs(distance) < abs(min_distance):
                min_distance = distance

    return distance


def bam_to_histogram(
        bam_path,
        add_to,
        feature_container,
        site_mode=False,
        bin_size=25,
        min_mq=30,
        max_distance=15_000,
        head=None,
        quick=True):
    histogram = add_to
    with pysam.AlignmentFile(bam_path) as f:
        i = 0
        for read in f:
            if head is not None and i > head:
                break
            if read.is_unmapped:
                continue
            if read.mapping_quality < min_mq:
                continue

            if read.has_tag('LY'):
                library = read.get_tag('LY')
            else:
                library = 'No_library_defined'
            chromosome = read.reference_name
            if quick or site_mode:
                if site_mode:
                    if not read.has_tag('DS'):
                        continue
                    i += 1
                    distance = distance_to_feature_start(chromosome, int(
                        read.get_tag('DS')), feature_container=feature_container)
                    if distance is not None and abs(distance) < max_distance:
                        histogram[library][np.floor(
                            distance / bin_size) * bin_size] += 1
                else:
                    i += 1
                    for distance in [
                        distance_to_feature_start(
                            chromosome,
                            read.reference_start,
                            feature_container=feature_container),
                        distance_to_feature_start(
                            chromosome,
                            read.reference_end,
                            feature_container=feature_container)]:

                        if distance is not None and abs(
                                distance) < max_distance:
                            histogram[library][np.floor(
                                distance / bin_size) * bin_size] += 1
            else:
                for q_pos, ref_pos in read.get_aligned_pairs(
                        matches_only=True, with_seq=False):
                    distance = distance_to_feature_start(
                        chromosome, ref_pos, feature_container=feature_container)
                    if distance is not None and abs(distance) < max_distance:
                        histogram[library][np.floor(distance / bin_size)] += 1
    print(f'{bam_path} read {i} reads')


if __name__ == '__main__':

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['figure.dpi'] = 160

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Visualise feature density of a bam file. (Coverage around stop codons, start codons, genes etc)')
    argparser.add_argument(
        '-o',
        type=str,
        help="output plot folder path, every library will be visualised separately",
        default='./plots/')
    argparser.add_argument(
        '-d',
        type=str,
        help="output data folder path, the data used for plotting will be exported to this folder",
        default='./tables/')
    argparser.add_argument(
        '-features',
        type=str,
        default='stop_codon,start_codon,exon,transcript',
        help="features to plot, separate by comma without space")
    argparser.add_argument('alignmentfiles', type=str, nargs='*')
    argparser.add_argument(
        '-gtf',
        type=str,
        required=True,
        help="GTF file containing the features to plot")
    argparser.add_argument(
        '-head',
        type=int,
        default=100_000_000,
        help="Use this amount of reads per bam file (or less if the bam file has less reads)")
    argparser.add_argument(
        '--bySite',
        action='store_true',
        help="Use the DS tag to count density")
    argparser.add_argument('--binSize', type=int, default=25)
    argparser.add_argument(
        '--maxDistance',
        type=int,
        default=15_000,
        help='Size of the window')
    argparser.add_argument(
        '--minMQ',
        type=int,
        default=30,
        help='Minimum mapping quality')
    args = argparser.parse_args()
    features = args.features.split(',')

    if not os.path.exists(args.o):
        os.makedirs(args.o)
    if not os.path.exists(args.d):
        os.makedirs(args.d)

    # Load the GTF files with the desired annotations
    annotations = {}
    for feature in features:

        annotations[feature] = singlecellmultiomics.features.FeatureContainer()
        annotations[feature].loadGTF(args.gtf, select_feature_type=[feature])
        if len(annotations[feature]) == 0:
            print(f"""Feature {feature} is not present in the suplied GTF file!
            Make sure to select feature types which are present in the GTF file!""")

    histograms = collections.defaultdict(
        lambda: collections.defaultdict(
            collections.Counter))  # feature->library->hist
    # Create histogram per feature / library
    for feature in features:
        if len(annotations[feature]) == 0:
            print(f"Skipping {feature} as it is not present in the GTF file")
            continue
        for bam_path in args.alignmentfiles:
            print(f"Now reading {bam_path} for annotation type {feature}")
            bam_to_histogram(
                bam_path,
                histograms[feature],
                annotations[feature],
                site_mode=args.bySite,
                bin_size=args.binSize,
                max_distance=args.maxDistance,
                min_mq=args.minMQ,
                head=args.head)

    # Write
    print('Writing tables')
    for feature, feature_data in histograms.items():
        for library, histogram in feature_data.items():
            df = pd.DataFrame(
                [list(histogram.keys()), list(histogram.values())]).transpose()
            df.columns = ['bin', 'observations']
            df.to_csv(
                f'{args.d}/{feature}_{library}_{args.maxDistance}_{args.binSize}.csv')

    # Plot the histograms
    print('Plotting')
    for feature, feature_data in histograms.items():
        for library, histogram in feature_data.items():
            print(f'{feature} {library}')
            fig, ax = plt.subplots(figsize=(12, 3))
            ax.scatter(
                list(
                    histogram.keys()), list(
                    histogram.values()), s=5, alpha=0.6)
            ax.set_xlabel(f'Distance from {feature} codon [bp]')
            ax.set_ylabel('# bases')
            ax.set_title(f'{library}')
            ax.set_ylim(bottom=-10)
            plt.tight_layout()
            plt.savefig(
                f'{args.o}/{feature}_{library}_{args.maxDistance}_{args.binSize}.png')
            plt.close()
