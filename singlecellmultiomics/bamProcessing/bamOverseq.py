#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
from collections import Counter, defaultdict
import pandas as pd
import seaborn as sns
from multiprocessing import Pool
from singlecellmultiomics.bamProcessing.bamBinCounts import blacklisted_binning_contigs
from more_itertools import grouper

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
            for i, read in enumerate(alignments.fetch(*location[:3])):
                if read.is_duplicate or read.mapping_quality < min_mq or read.is_read2 or read.is_qcfail or (allelic and not read.has_tag(
                        'DA')):
                    continue

                bin_index = int(read.get_tag('DS') / bin_size)

                location_key = (
                        read.reference_name,
                        bin_index * bin_size,
                        (bin_index + 1) * bin_size
                    )

                if allelic:
                    location_key = (read.get_tag('DA'), *location_key)

                overseq[location_key ][read.get_tag('SM')][read.get_tag('af')] += 1
                # if i%1_000_000==0:
                #    print(i,read.reference_name, read.reference_start,end='\r')

    return overseq

def obtain_overseq_dictionary(bam_path, bin_size, min_mq=10, allelic=False):
    overseq = defaultdict(overseq_dict)

    worker_bins = list(grouper(blacklisted_binning_contigs(contig_length_resource=bam_path,
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

