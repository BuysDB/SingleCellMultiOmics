#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
from collections import Counter, defaultdict
import pandas as pd
import seaborn as sns
from multiprocessing import Pool
from singlecellmultiomics.bamProcessing.bamBinCounts import blacklisted_binning_contigs
from more_itertools import grouper
import argparse
import os

def overseq_dict():
    return defaultdict(Counter)

def merge_overseq_dicts(into, source):
    for k, overseq_per_sample in source.items():
        for sample, overseq_for_sample in overseq_per_sample.items():
            into[k][sample] += overseq_for_sample

def create_overseq_distribution(args):
    locations, bin_size, bam_path, min_mq, allelic = args
    overseq = defaultdict(overseq_dict)
    with pysam.AlignmentFile(bam_path) as alignments:
        for location in locations:
            if location is None:
                continue

            # Location is a tuple (contig, start, end, fetch_start, fetch_end)
            # We only need (contig, start, end) to obtain reads in the target region
            for i, read in enumerate(alignments.fetch(*location[:3])):

                start, end = location[-2:]
                if read.is_duplicate or read.mapping_quality < min_mq or read.is_read2 or read.is_qcfail or (allelic and not read.has_tag(
                        'DA')):
                    continue

                # Prevent double counting of fragments
                if read.has_tag('DS'):
                    site = read.get_tag('DS')
                    if site < start or site > end:
                        continue
                    # @todo ; implement for reads without DS tag

                bin_index = int(read.get_tag('DS') / bin_size)

                location_key = (
                        read.reference_name,
                        bin_index * bin_size,
                        (bin_index + 1) * bin_size
                    )

                if allelic:
                    # Prepend allele to bin identifier:
                    location_key = (read.get_tag('DA'), *location_key)

                overseq[location_key ][read.get_tag('SM')][read.get_tag('af')] += 1
                # if i%1_000_000==0:
                #    print(i,read.reference_name, read.reference_start,end='\r')

    return overseq

def obtain_overseq_dictionary(bam_path, bin_size, min_mq=10, allelic=False, blacklist_path=None):
    """
    Create molecule duplicate counting dictionary

    Args:
        bam_path: path to bam file to obtain duplicate counts from (requires af tag to be set)
        bin_size: bin size for output dict
        min_mq: minimum mapping quality
        allelic: wether to split bins into multiple alleles based on the DA tag value
        blacklist_path: path to blacklist bed file (optional)

    Returns:
        overseq (dict) : {(location):{sample(str):{reads_per_molecule(int):count (int)}}}

    """



    overseq = defaultdict(overseq_dict)

    worker_bins = list(grouper(blacklisted_binning_contigs(contig_length_resource=bam_path,
                                                           blacklist_path=blacklist_path,
                                                           bin_size=bin_size,
                                                           fragment_size=0), 50))  # take N bins for every worker

    with Pool() as workers:
        for i, overseq_for_bin in enumerate(workers.imap_unordered(create_overseq_distribution, (
                (
                        locations,
                        bin_size,
                        bam_path,
                        min_mq,
                        allelic
                )

                for locations in worker_bins))):

            merge_overseq_dicts(overseq, overseq_for_bin)

            print(f'{((i / len(worker_bins)) * 100):.2f} % ({i}/{len(worker_bins)})       ', end='\r')

    return overseq

def write_overseq_dict_to_single_sample_files(overseq, target_dir):
    # reads per cell:
    reads_per_cell = Counter()
    locations = sorted(list(overseq.keys()))

    for k, overseq_per_cell in overseq.items():
        for s in overseq_per_cell:
            reads_per_cell[s] += sum([copies * seen for copies, seen in overseq_per_cell[s].items()])
    reads_per_cell_dict = reads_per_cell
    reads_per_cell = pd.DataFrame({'reads': reads_per_cell})
    selected_cells = reads_per_cell[reads_per_cell['reads'] > 100].index
    selected_cells = sorted(selected_cells)

    for sample in selected_cells:
        cell_hist = {}
        cell_map = defaultdict(dict)
        for ki, k in enumerate(locations):
            overseq_per_cell = overseq[k]
            overseq_counter_for_cell = overseq_per_cell[sample]
            cell_hist[k] = overseq_counter_for_cell

        print(sample)
        pd.DataFrame(cell_hist).sort_index().sort_index(1).to_csv(f'{target_dir}/cell_hist_{sample}.csv')




if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Create oversequencing tables for single cells""")
    argparser.add_argument('bamfile', metavar='bamfile', type=str, help='Bam file where duplicates have been flagged \
    , and the amount of reads associated to the duplicate stored in the sam tag "af"')
    argparser.add_argument('-output_folder', default='./overseq_tables', help='Folder to write sample CSV files to')
    argparser.add_argument('--allelic', action='store_true', help='Split bins into alleles (based on DA tag)')
    argparser.add_argument('-bin_size', type=int, default=500_000, help='Bin size for the output files')
    argparser.add_argument('-blacklist', type=str, help='Bed file containing regions to ignore from counting')

    fi = argparser.add_argument_group("Filters")
    fi.add_argument('-min_mapping_qual', default=40, type=int)

    args = argparser.parse_args()

    # Create output folder if it does not exist:
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    overseq = obtain_overseq_dictionary(args.bamfile, args.bin_size, min_mq=args.min_mapping_qual, allelic=args.allelic, blacklist_path=args.blacklist)
    # (location) : cell : duplicates : count
    write_overseq_dict_to_single_sample_files(overseq, args.output_folder)
