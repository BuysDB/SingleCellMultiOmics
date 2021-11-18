#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pysam
import numpy as np
import os
import pickle
import gzip
import pandas as pd
import multiprocessing
from singlecellmultiomics.bamProcessing import get_contig_sizes, get_contig_size, get_contigs
from statsmodels.nonparametric.smoothers_lowess import lowess
from datetime import datetime
from itertools import chain
from more_itertools import windowed
from typing import Generator
from singlecellmultiomics.methylation import MethylationCountMatrix
from pysamiterators import CachedFasta
from pysam import FastaFile
from singlecellmultiomics.utils import reverse_complement, is_main_chromosome, pool_wrapper
from collections import defaultdict, Counter
from itertools import product
from singlecellmultiomics.bamProcessing.bamFunctions import mate_iter
from multiprocessing import Pool
from singlecellmultiomics.pyutils import sorted_slice


def _generate_count_dict(args):
    """
    Obtain counts for a bam path, given a bin size and region

    args:
        args (tuple) : bam_path, bin_size, contig, start, stop
    """
    bam_path, bin_size, contig, start, stop, filter_function = args #reference_path  = args

    #reference_handle = pysam.FastaFile(reference_path)
    #reference = CachedFasta(reference_handle)

    cut_counts = defaultdict(Counter )
    i = 0
    with pysam.AlignmentFile(bam_path) as alignments:

        for R1,R2 in mate_iter(alignments, contig=contig, start=start, stop=stop):

            if filter_function is None:
                if R1 is None or R1.is_duplicate or R1.is_qcfail:
                    continue
            else:
                if not filter_function(R1,R2):
                    continue

            if not R1.has_tag('DS'):
                cut_pos = R1.reference_end if  R1.is_reverse else R1.reference_start
            else:
                cut_pos = int(R1.get_tag('DS'))

            if cut_pos is None:
                continue

            if R1.has_tag('SM'):
                sample = R1.get_tag('SM')
            else:
                sample = 'No_Sample'

            bin_idx=int(cut_pos/bin_size)*bin_size
            cut_counts[(contig,bin_idx)][sample] += 1

    return cut_counts, contig, bam_path

# cell-> context -> obs
def _generate_count_dict_prefixed(args):

    bam_path, bin_size, contig, start, stop, filter_function,  alias, prefix = args #reference_path  = args

    cut_counts = defaultdict(Counter )
    i = 0
    with pysam.AlignmentFile(bam_path) as alignments:

        for R1,R2 in mate_iter(alignments, contig=contig, start=start, stop=stop):

            if filter_function is None:
                if R1 is None or R1.is_duplicate or not R1.has_tag('DS') or R1.is_qcfail:
                    continue
            else:
                if not filter_function(R1,R2):
                    continue

            cut_pos = R1.get_tag('DS')
            sample = R1.get_tag('SM')

            if prefix is not None:
                sample = (prefix, sample)

            if alias is not None:

                cut_counts[alias][sample] += 1
            else:
                bin_idx=int(cut_pos/bin_size)*bin_size

                cut_counts[(contig,bin_idx)][sample] += 1

    return cut_counts, contig, bam_path

# Get TA fraction per cell
def _generate_ta_count_dict_prefixed(args):

    bam_path, bin_size, contig, start, stop, filter_function, alias, prefix = args

    cut_counts = defaultdict(Counter )
    i = 0
    with pysam.AlignmentFile(bam_path) as alignments:

        for R1,R2 in mate_iter(alignments, contig=contig, start=start, stop=stop):

            if R1 is None or R1.is_duplicate or not R1.has_tag('DS') or R1.is_qcfail:
                continue
            if R1.get_tag('lh')!='TA':
                continue

            if filter_function is not None and not filter_function(R1,R2):
                continue

            cut_pos = R1.get_tag('DS')
            sample = R1.get_tag('SM')
            if prefix is not None:
                sample = (prefix, sample)

            if alias is not None:

                cut_counts[alias][sample] += 1
            else:
                bin_idx=int(cut_pos/bin_size)*bin_size

                cut_counts[(contig,bin_idx)][sample] += 1

    return cut_counts, contig, bam_path



def _get_contig_cuts(bam: str, contig: str, mqthreshold=50):
    """
    Obtain np.array of all cuts for the supplied contig and bam file for each sample indiviually,
    the cut arrays are sorted by genomic coordinate

    Returns:
        contig, {cell: np.array([cut1(int),cut2(int)])}
    """
    single_cell_cuts = defaultdict(list)
    with pysam.AlignmentFile(bam) as a:
        for read in a.fetch(contig):
            if read.is_read1 and not read.is_qcfail and not read.is_duplicate and read.mapping_quality>mqthreshold:
                single_cell_cuts[read.get_tag('SM')].append(read.get_tag('DS'))

    return contig, {

        cell:np.array(sorted(cuts))
        for cell, cuts in single_cell_cuts.items()
    }


def bam_to_single_sample_cuts(bams: list, n_processes=None, contigs=None, mqthreshold=50):
    """
    Obtain np.array of all cuts for the supplied contig and bam file for each sample indiviually,
    the cut arrays are sorted by genomic coordinate
    When contigs is not supplied cuts are determined for all contigs except random scaffolds and alternative alleles

    Returns:
        cell_cuts: {cell: {contig: np.array([cut1(int),cut2(int)])}},
        total_cuts_per_cell : {cell:total_cuts (int)}
    """

    cell_cuts = {} # cell -> contig -> cuts
    total_cuts_per_cell = Counter()
    with Pool(n_processes) as workers:
        if contigs is None:
            contigs = [c for c in get_contigs(bams[0]) if is_main_chromosome(c)]


        for contig, cuts in workers.imap(pool_wrapper, ((
            _get_contig_cuts, {'contig':contig,
                               'bam':bam,
                               'mqthreshold':mqthreshold
                               })

            for contig, bam in  product(contigs,bams))):

            for cell, cc in cuts.items():
                if not cell in cell_cuts:
                    cell_cuts[cell] = {}
                cell_cuts[cell][contig] = cc
                total_cuts_per_cell[cell] += len(cc)

    return cell_cuts, total_cuts_per_cell


def get_binned_counts_prefixed(bam_dict, bin_size, regions=None, ta=False, filter_function=None, n_threads=None) -> pd.DataFrame:
    """
    Same as get_binned_counts but accespts a bam_dict, {'alias':[bam file, bam file], 'alias2':[bam file, ]}
    The alias is added as the first level of the output dataframe multi-index

    Returns:
        pd.DataFrame
    """
    fs = 1000

    aliased = False

    if regions is None:
        regions = [(c,None,None) for c in get_contig_sizes(list(bam_dict.values())[0][0]).keys()]

    else:
        for i,r in enumerate(regions):
            if type(r)==str:
                regions[i] = (r,None,None)

            elif len(r)==3:

                contig, start, end =r
                if type(start)==int:
                    start = max(0,start-fs)

                regions[i] = (contig,start,end)


            else:
                contig, start, end, alias =r
                if type(start)==int:
                    start = max(0,start-fs)

                regions[i] = (contig,start,end, alias)
                aliased=True

    #if aliased:
    #    jobs = [(bam_path, bin_size, *region) for region, bam_path in product(regions, bams)]

    jobs = []
    for bam_group, bam_list in bam_dict.items():
        for region in regions:
            for bam_path in bam_list:
                jobs.append( (bam_path, bin_size, *region, filter_function, None, bam_group) )
            #jobs = [(bam_path, bin_size, *region, None, prefix) for region, bam_path in product(regions, bams)]


    cut_counts = defaultdict(Counter)
    with Pool(n_threads) as workers:

        for i, (cc, contig, bam_path) in enumerate(workers.imap(_generate_count_dict_prefixed if not ta else _generate_ta_count_dict_prefixed,jobs)):

            for k,v in cc.items():
                cut_counts[k] += v

            print(i,'/', len(jobs), end='\r')

    return pd.DataFrame(cut_counts).T


def get_binned_counts(bams, bin_size, regions=None, filter_function=None, n_threads=None):

    fs = 1000
    if regions is None:
        regions = [(c,None,None) for c in get_contig_sizes(bams[0]).keys()]

    else:
        for i,r in enumerate(regions):
            if type(r)==str:
                regions[i] = (r,None,None)
            else:
                contig, start, end =r
                if type(start)==int:
                    start = max(0,start-fs)

                regions[i] = (contig,start,end)

    jobs = [(bam_path, bin_size, *region, filter_function) for region, bam_path in product(regions, bams)]


    cut_counts = defaultdict(Counter)
    with multiprocessing.Pool(n_threads) as workers:

        for i, (cc, contig, bam_path) in enumerate(workers.imap(_generate_count_dict,jobs)):

            for k,v in cc.items():
                cut_counts[k] += v

            print(i,'/', len(jobs), end='\r')

    return pd.DataFrame(cut_counts).T


def fill_range(start, end, step):
    """
    range iterator from start to end with stepsize step
    but generates a final step which steps to end

    Args:
        start (int) : start

        end (int)

        step (int) : step


    Example:
        >>> list(fill_range(0,10,4))
        [(0, 4), (4, 8), (8, 10)]


    """

    e = start
    for s in range(start, end, step):
        e = s + step
        if e > end:
            e = e - step
            break
        yield s, e

    if e < end:
        yield e, end


def trim_rangelist(rangelist, start, end):
    """
    Trim list of start, end coordinates to only keep ranges within start and end, and trim ranges which overlap with start/end.

    Args:

        rangelist(list) : list of tuples (start,end)

        start(int) : inclusive start coordiante

        end(int) : exclusive end coordinate

    Yields:
        (range_start, range_end)

    """
    for s, e in rangelist:

        overlap = False
        if s >= start and s < end:
            overlap = True
        if e >= start and e < end:
            overlap = True

        if not overlap:
            continue

        yield max(s, start), min(e, end)


def get_bins_from_bed_iter(path, contig=None):

    with gzip.open(path) if path.endswith('.gz') else open(path)  as f:
        for line in f:
            c, start, end = line.strip().split(None, 3)[:3]
            start, end = int(start), int(end)
            if contig is None or c == contig:
                yield c, start, end


def get_bins_from_bed_dict(path, contig=None):
    bd = {}
    for c, start, end in get_bins_from_bed_iter(path, contig=contig):
        if not c in bd:
            bd[c] = list()
        bd[c].append((start, end))

    return bd


def range_contains_overlap(clist):
    """
    Check if the supplied iteratorof start,end coordinates contains an overlap

    Args:
        clist (list) : sorted list of start,end coordinates of regions
    """
    clist = sorted(clist)
    # check for overlaps...
    if len(clist) < 2:
        return False

    for (start, end), (next_start, next_end) in windowed(clist, 2):
        if start > next_start or end > next_start or end > next_end or start > next_end:
            # print(f'Overlap between {start}:{end} and {next_start}:{next_end}')
            return True
    return False


def _merge_overlapping_ranges(clist):
    merged = False
    for (start, end), (next_start, next_end) in windowed(clist, 2):
        if merged:
            merged = False
            continue

        if start > next_start or end > next_start or end > next_end or start > next_end:
            yield (min(start, next_start), max(next_end, end))
            merged = True

        else:
            yield start, end
    if not merged:
        yield clist[-1]


def merge_overlapping_ranges(clist):
    clist = sorted(clist)
    while range_contains_overlap(clist):
        clist = sorted(list(_merge_overlapping_ranges(clist)))
    return clist


def blacklisted_binning_contigs(contig_length_resource: str, bin_size: int, fragment_size: int,
                                blacklist_path: str = None, contig_whitelist: list = None) -> Generator:
    """
    Generate a list of (contig, bin_start, bin_end) tuples of size bin_size or smaller or when fragment_size is supplied
    (contig, bin_start, bin_end, fetch_start, fetch_end). All regions present in the blacklist BED file will not be
    part of the generated bins.

    Args:
        contig_length_resource(str): Path to bam file from which to extract the contig lengths
        bin_size(int) : maximum size of generated bins (might produce some bins which are smaller)
        fragment_size(int) : When this value is supplied fetch_start, fetch_end will be produced which will be equal to
                            bin_start-fragment_size and bin_end+fragment size. But will never overlap with blacklisted
                            regions or exceed contig boundaries.
        blacklist_path(str): path to blacklist bed file
        contig_whitelist(iterable): A set of contigs to only include in the result. All contigs are included when
                                    contig_whitelist is not specified.

    Returns:
        bin_tuples(Generator): (contig, bin_start, bin_end),
                               ( contig, bin_start, bin_end, fetch_start, fetch_end ) when fragment_size is specified
    """

    if blacklist_path is not None:
        blacklist_dict = get_bins_from_bed_dict(blacklist_path)
    else:
        blacklist_dict = {}

    for contig, length in (get_contig_sizes(contig_length_resource).items()
    if type(contig_length_resource) is str else contig_length_resource):
        if contig_whitelist is not None and not contig in contig_whitelist:
            continue

        if fragment_size is not None:
            for bin_start, bin_end, fetch_start, fetch_end in \
                    blacklisted_binning(
                        start_coord=0,
                        end_coord=length,
                        bin_size=bin_size,
                        blacklist=sorted(blacklist_dict.get(contig, [])),
                        fragment_size=fragment_size):
                yield contig, bin_start, bin_end, fetch_start, fetch_end
        else:
            for bin_start, bin_end in \
                    blacklisted_binning(
                        start_coord=0,
                        end_coord=length,
                        bin_size=bin_size,
                        blacklist=sorted(blacklist_dict.get(contig, []))):
                yield contig, bin_start, bin_end


def blacklisted_binning(start_coord: int, end_coord: int, bin_size: int, blacklist: list = None,
                        fragment_size: int = None):
    """
    Obtain a list of regions to fetch between start coord and end_coord with given bin_size and
    excluding the regions in the blacklist. Optimizes the bin size and allows for a fragment_size parameter to be set
    which expands the bins with fragment size without overlap with blacklist regions and
    start_coord, end_coord boundaries

    Args:
        start_coord(int)

        end_coord(int)

        bin_size(int)

        blacklist(list) : format (start,end), (start,end)

        fragment_size(int)

    Yields:
        bin_start, bin_end when fragment_size is None, bin_start, bin_end, fetch_start, fetch_end otherwise

    """
    if blacklist is None:
        blacklist = []
    elif len(blacklist) > 1:
        blacklist = merge_overlapping_ranges(blacklist)
    current = start_coord
    for i, (start, end) in enumerate(chain(
            trim_rangelist(blacklist, start_coord, end_coord),
            [(end_coord, end_coord + 1), ])):

        # Fill from current to start:
        if start == current:
            current = end
            continue

        total_bins = len(list(fill_range(current, start, bin_size)))
        if total_bins < 0:
            continue
        if total_bins == 0:
            total_bins = 1

        local_bin_size = int((start - current) / total_bins)

        for fi, (pos_s, pos_e) in enumerate(fill_range(current, start, local_bin_size)):

            if fragment_size is None:
                yield pos_s, pos_e
            else:
                fs = pos_s - fragment_size
                fe = pos_e + fragment_size
                if fi == 0:
                    # Left clip, no available bins on the left
                    fs = pos_s

                if fi == total_bins - 1:
                    fe = pos_e

                yield pos_s, pos_e, fs, fe
            current = pos_e
        current = end


def obtain_approximate_reference_cut_position(site: int, contig: str, alt_spans: dict) -> tuple:
    alt_contig, alt_start, alt_end = alt_spans[contig]
    return contig, site + alt_start


def obtain_counts(commands, reference, live_update=True, show_n_cells=4, update_interval=3, threads=4, count_function=None, show_progress=False):
    if count_function is None:
        count_function = count_fragments_binned

    if live_update:
        from singlecellmultiomics.utils.plotting import GenomicPlot
        import matplotlib.pyplot as plt
        cell_plots = {}
        for cell_index in range(show_n_cells):
            gplot = GenomicPlot(reference)
            fig = gplot.get_figure()
            fig.canvas.draw()
            cell_plots[cell_index] = {'plot': gplot, 'fig': fig}

        plt.pause(0.01)

    # import random
    # random.shuffle(commands)

    counts = {}

    prev = None

    top_cells = None

    start_time = datetime.now()

    update_method = 'partial_df'

    if show_progress:
        commands = list(commands) # generator -> list
        print(f"Counting started ({len(commands)})", end = '\r')

    with multiprocessing.Pool(threads) as workers:

        for i, result in enumerate(workers.imap_unordered(count_function,
                                                          commands)):

            if show_progress:
                print(f"Counting progress { ((i/len(commands))*100):.2f} % ({i} / {len(commands)}) ", end = '\r')
            # Result is of the format: counts[bin_id][sample] = obs (int)
            #counts.update(result)
            for bin_id, sample_dict in result.items():
                if not bin_id in counts:
                    counts[bin_id] =  sample_dict
                else:
                    counts[bin_id].update(sample_dict)



            if live_update and update_method == 'partial_df':
                if (datetime.now() - start_time).total_seconds() > 2 and (
                        prev is None or (datetime.now() - prev).total_seconds() >= update_interval):
                    if len(result) == 0:
                        continue

                    df = pd.DataFrame(counts).T
                    if df.sum().sum() == 0:
                        continue
                    prev = datetime.now()

                    if top_cells is None:
                        top_cells = df.sum().sort_values()[-show_n_cells:]

                    df = df[top_cells.index].fillna(0)
                    df = np.clip(0, 2, df / np.percentile(df, 99, axis=0))

                    for contig in [list(result.keys())[0][0]]:
                        x = np.array([(stop + start) / 2 for start, stop in df.loc[contig].index.values])

                        for cell_index, (cell, row) in enumerate(df.loc[contig].T.iterrows()):
                            gplot = cell_plots[cell_index]['plot']
                            gplot.reset_axis(contig)
                            fig = cell_plots[cell_index]['fig']
                            fig.suptitle(cell)
                            gplot[contig].scatter(x, row.values, s=0.1, c='k')
                            fig.canvas.draw()
                    plt.pause(0.001)

    if show_progress:
        print("Counting finished     ")

    # Show final result
    if live_update:
        df = pd.DataFrame(counts).T
        df = df[top_cells.index].fillna(0)
        df = np.clip(0, 2, df / np.percentile(df, 99, axis=0))

        for contig in cell_plots[0]['plot'].contigs:
            x = np.array([(stop + start) / 2 for start, stop in df.loc[contig].index.values])

            for cell_index, (cell, row) in enumerate(df.loc[contig].T.iterrows()):
                gplot = cell_plots[cell_index]['plot']
                gplot.reset_axis(contig)
                fig = cell_plots[cell_index]['fig']
                fig.suptitle(cell)
                gplot[contig].scatter(x, row.values, s=0.1, c='k')
                fig.canvas.draw()
        plt.pause(0.001)

    return counts


def read_counts(read, min_mq, dedup=True, read1_only=False,ignore_mp=False, ignore_qcfail=False, verbose=False):

    if read1_only and (not read.is_read1 or read is None):
        if verbose:
            print('NOT READ1')
        return False

    if read.is_qcfail and not ignore_qcfail:
        if verbose:
            print('QCFAIL')
        return False

    if dedup and read.is_duplicate:
        if verbose:
            print('DUPLICATE')
        return False
    if not ignore_mp and read.has_tag('mp') and read.get_tag('mp') != 'unique':
        if verbose:
            print('MP')
        return False
    if (min_mq is not None and read.mapping_quality < min_mq):
        if verbose:
            print('MAPQ')
        return False

    if verbose:
        print('OK')
    # if read.has_tag('RZ') and read.get_tag('RZ') != 'CATG':
    #    return False

    return True


def gc_correct(args):
    observations, gc_vector, MAXCP = args
    correction = lowess(observations, gc_vector)
    return np.clip(observations / np.interp(gc_vector, correction[:, 0], correction[:, 1]), 0, MAXCP)


def genomic_bins_to_gc(reference, df):
    bins_to_gc = {}
    for k in df.columns:
        if len(k)==4:
            (allele, contig, start, end) = k
        else:
            (contig, start, end) = k
        if not k in bins_to_gc:
            sequence = reference.fetch(contig, start, end).upper()
            gc = sequence.count('G') + sequence.count('C')
            div = ((sequence.count('A') + sequence.count('T') + gc))
            if div == 0:
                # There is no data, plop the mean gc in later
                bins_to_gc[k] = np.nan
            else:
                gcrat = (gc) / div
                bins_to_gc[k] = gcrat
    return bins_to_gc


def gc_correct_cn_frame(df, reference, MAXCP, threads, norm_method='median'):
    # Perform GC correction
    chrom_sizes = dict(zip(reference.references, reference.lengths))
    # Extract GC percentage from reference for the selected bin size:

    bins_to_gc =  genomic_bins_to_gc(reference, df)
    qf = pd.DataFrame({'gc': bins_to_gc})
    qf = qf.fillna(qf['gc'].mean())
    # Join the GC table with the count matrix
    gc_matched = df.T.join(qf, how='left')['gc']

    # This performs GC correction for every cell using loess regression
    with multiprocessing.Pool(threads) as workers:
        keep_bins = df.columns
        gc_vector = gc_matched[keep_bins]

        corrected_cells = list(workers.imap(
            gc_correct, [(row, gc_vector.values, MAXCP) for cell, row in df.iterrows()]))

    corrected_cells = pd.concat(corrected_cells, axis=1).T
    if norm_method == 'median':
        corrected_cells = ((corrected_cells.T / corrected_cells.median(1)) * 2).T
    elif norm_method == 'mean':
        corrected_cells = ((corrected_cells.T / corrected_cells.mean(1)) * 2).T
    else:
        raise ValueError('norm_method not understood')

    return corrected_cells


def generate_jobs(alignments_path, bin_size=1_000_000, bins_per_job=10):
    for job_group in (((contig, start, start + bin_size * bins_per_job)
                       for start in range(0, length, bin_size * bins_per_job))
                      for contig, length in
                      get_contig_sizes(alignments_path).items()):
        yield from job_group


def generate_commands(alignments_path,
                      bin_size=1_000_000,
                      bins_per_job=10,
                      alt_spans=None,
                      min_mq=50,
                      max_fragment_size=1000,
                      head=None,
                      key_tags=None,
                      dedup=True,
                      kwargs=None,
                      skip_contigs=None
                      ):


    if type(alignments_path) is list:
        iterfiles =  alignments_path
    else:
        iterfiles = [alignments_path]

    for alignments_path in iterfiles:
        for i, (contig, start, end) in enumerate(
                generate_jobs(alignments_path=alignments_path, bin_size=bin_size, bins_per_job=bins_per_job)):

            if skip_contigs is not None and contig in skip_contigs:
                continue
            yield (alignments_path, bin_size, max_fragment_size, \
                   contig, start, end, \
                   min_mq, alt_spans, key_tags, dedup, kwargs)
            if head is not None and i >= (head - 1):
                break


def count_methylation_binned(args):
    (alignments_path, bin_size, max_fragment_size, \
     contig, start, end, \
     min_mq, alt_spans, key_tags, dedup, kwargs) = args


    # single_location => count single cpgs
    single_location = kwargs.get('single_location', False)
    dyad_mode =  kwargs.get('dyad_mode', False)

    # Stranded mode:
    stranded = kwargs.get('stranded', False)

    contexts_to_capture = kwargs.get('contexts_to_capture', None) # List

    count_reads = kwargs.get('count_reads', False)

    reference = None

    if contexts_to_capture is not None:
        context_radius = kwargs.get('context_radius', None)
        reference_path = kwargs.get('reference_path', None)
        if reference_path is not None:
            reference = CachedFasta(FastaFile(reference_path))



    default_sample_name = kwargs.get('default_sample_name', 'bulk_sample')

    min_counts_per_bin = kwargs.get('min_counts_per_bin',10) # Min measurements across all cells
    # Cant use defaultdict because of pickles :\
    met_counts = MethylationCountMatrix(threads=kwargs.get('threads', None))  # Sample->(contig,bin_start,bin_end)-> [methylated_counts, unmethylated]

    if count_reads:
        read_count_dict = defaultdict(Counter)  # location > sample > read_obs

    # Define which reads we want to count:
    known =  set()
    if 'known' in kwargs and kwargs['known'] is not None:
        # Only ban the very specific TAPS conversions:
        try:
            with pysam.VariantFile(kwargs['known']) as variants:
                for record in variants.fetch(contig, start, end):
                    if record.ref=='C' and 'T' in record.alts:
                        known.add( record.pos)
                    if record.ref=='G' and 'A' in record.alts:
                        known.add(record.pos)
        except ValueError:
            # This happends on contigs not present in the vcf
            pass

    p = 0
    start_time = datetime.now()
    with pysam.AlignmentFile(alignments_path, threads=4) as alignments:
        # Obtain size of selected contig:
        contig_size = get_contig_size(alignments, contig)
        if contig_size is None:
            raise ValueError('Unknown contig')

        # Determine where we start looking for fragments:
        f_start = max(0, start - max_fragment_size)
        f_end = min(end + max_fragment_size, contig_size)

        for p, read in enumerate(alignments.fetch(contig=contig, start=f_start,
                                                  stop=f_end)):

            if p%50==0 and 'maxtime' in kwargs and kwargs['maxtime'] is not None:
                if (datetime.now() - start_time).total_seconds() > kwargs['maxtime']:
                    print(f'Gave up on {contig}:{start}-{end}')
                    break

            if not read_counts(read, min_mq=min_mq, dedup=dedup, verbose=False):
                continue

            tags = dict(read.tags)
            sample = tags.get('SM', default_sample_name)

            if count_reads and (read.is_read1 or not read.is_paired):

                site = tags.get('DS', (read.reference_end if read.is_reverse else read.reference_start))
                # Obtain the bin index
                if single_location:
                    bin_start = site
                    bin_end = site + 1
                else:
                    bin_i = int(site / bin_size)
                    bin_start = bin_size * bin_i
                    bin_end = min(bin_size * (bin_i + 1), contig_size)

                read_count_dict[(read.reference_name,bin_start, bin_end)][sample] += 1

            for i, (qpos, site) in enumerate(read.get_aligned_pairs(matches_only=True)):

                # Don't count sites outside the selected bounds
                if site < start or site >= end:
                    continue

                if site in known:
                    continue

                call = tags['XM'][i]

                final_call = None
                if contexts_to_capture is not None:
                    if call=='.':
                        continue

                    try:
                        context = reference.fetch(read.reference_name,
                                                  site - context_radius,
                                                  site + context_radius + 1).upper()
                    except KeyError:
                        raise KeyError(f'The supplied reference file does not match the reads. Contig "{read.reference_name}" missing from supplied reference.')

                    if len(context)!=(context_radius*2 +1):
                        continue

                    if context[context_radius] == 'G':
                        context = reverse_complement(context)

                    if context in contexts_to_capture:
                        final_call = call.isupper()

                    else:
                        continue

                    if context[context_radius] not in 'CG':
                        # This is bad... the position is not  a C or G. It might be another IUPAC nucleotide code.

                        #Ditch this READ if the ref base is A/T
                        if context[context_radius] in 'AT':
                            print(f'For read {read.reference_name}:{read.cigarstring} at {site}, query pos {qpos} a {context[context_radius]} was found in the reference. Ignoring this read')
                            break

                        #Ditch this BASE if IUPAC nucleotide code:
                        continue


                elif call in 'Zz':

                    final_call = call=='Z'
                    #CpG
                    #^
                    if read.is_reverse and dyad_mode:
                        #GpC
                        #  ^
                        site += 1

                else:
                    continue


                # Process alternative contig counts:
                if alt_spans is not None and contig in alt_spans:
                    contig, site = obtain_approximate_reference_cut_position(
                        site, contig, alt_spans)

                # Obtain the bin index
                if single_location:
                    bin_i = site
                    bin_start = site
                    bin_end = site + 1
                else:
                    bin_i = int(site / bin_size)
                    bin_start = bin_size * bin_i
                    bin_end = min(bin_size * (bin_i + 1), contig_size)

                # Add additional tag information: (For example the allele tag)
                if key_tags is not None:
                    tag_values = [(read.get_tag(tag) if read.has_tag(tag) else None) for tag in key_tags]
                    if stranded:
                        bin_id = (*tag_values, contig, bin_start, bin_end, '+-'[read.is_reverse])
                    else:
                        bin_id = (*tag_values, contig, bin_start, bin_end)
                else:

                    if stranded:
                        bin_id = (contig, bin_start, bin_end,  '+-'[read.is_reverse])
                    else:
                        bin_id = (contig, bin_start, bin_end)



                if final_call is not None:
                    met_counts[sample, bin_id][final_call]+=1

    met_counts.prune(min_samples=kwargs.get('min_samples',0), min_variance=kwargs.get('min_variance',0))
    if reference is not None:
        reference.handle.close()

    return met_counts, (read_count_dict if count_reads else None)



def count_fragments_binned(args):
    (alignments_path, bin_size, max_fragment_size, \
     contig, start, end, \
     min_mq, alt_spans, key_tags, dedup, kwargs) = args

    counts = {}  # Sample->(contig,bin_start,bin_end)->counts
    # Define which reads we want to count:

    p = 0
    ignore_mp = kwargs.get('ignore_mp',False)
    with pysam.AlignmentFile(alignments_path, threads=4) as alignments:
        # Obtain size of selected contig:
        contig_size = get_contig_size(alignments, contig)
        if contig_size is None:
            raise ValueError('Unknown contig')

        # Determine where we start looking for fragments:
        f_start = max(0, start - max_fragment_size)
        f_end = min(end + max_fragment_size, contig_size)

        for p, read in enumerate(alignments.fetch(contig=contig, start=f_start,
                                                  stop=f_end)):

            if not read_counts(read, min_mq=min_mq, dedup=dedup, read1_only=True,ignore_mp=ignore_mp):
                continue

            # Extract the site
            site = int(read.get_tag('DS'))

            # Don't count sites outside the selected bounds
            if site < start or site >= end:
                continue

            sample = read.get_tag('SM')

            # Process alternative contig counts:
            if alt_spans is not None and contig in alt_spans:
                contig, site = obtain_approximate_reference_cut_position(
                    site, contig, alt_spans)

            # Obtain the bin index
            bin_i = int(site / bin_size)
            bin_start = bin_size * bin_i
            bin_end = min(bin_size * (bin_i + 1), contig_size)

            # Add additional tag information: (For example the allele tag)
            if key_tags is not None:
                tag_values = [(read.get_tag(tag) if read.has_tag(tag) else None) for tag in key_tags]
                bin_id = (*tag_values, contig, bin_start, bin_end)
            else:
                bin_id = (contig, bin_start, bin_end)

            # Add a (single) count tot the dictionary:
            if not bin_id in counts:
                counts[bin_id] = {}

            if not sample in counts[bin_id]:
                counts[bin_id][sample] = 1
            else:
                counts[bin_id][sample] += 1

    return counts


def count_fragments_binned_wrap(args):
    (alignments_path, bin_size, max_fragment_size, \
     contig, start, end, \
     min_mq, alt_spans, dedup, kwargs) = args

    tp = f'./TEMP_{contig}_{start}.pickle.gz'
    res = os.system(f'bamBinCounts.py {alignments_path} -o {tp} -start {start} -end {end} -contig {contig}')
    with gzip.open(tp, 'rb') as pf:
        result = pickle.load(pf)
    os.remove(tp)
    return result


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract counts from bam file')

    argparser.add_argument('alignmentfile', type=str)
    argparser.add_argument('-head', type=int, default=None)
    argparser.add_argument('-bin_size', type=int, default=250_000)
    argparser.add_argument('-t', type=int, default=16, help='Threads')
    argparser.add_argument('-j', type=int, default=5, help='Bins per worker')
    argparser.add_argument(
        '-min_mq',
        type=int,
        default=30,
        help='Minimum mapping quality')
    argparser.add_argument('-max_fragment_size', type=int, default=1000)
    argparser.add_argument('-o', type=str, required=True)
    args = argparser.parse_args()

    counts = {}

    with multiprocessing.Pool(args.t) as workers:

        job_gen = generate_commands(
            alignments_path=args.alignmentfile,
            bin_size=args.bin_size,
            bins_per_job=args.j,
            head=args.head)

        for job_total, _ in enumerate(job_gen):
            pass

        job_gen = generate_commands(
            alignments_path=args.alignmentfile,
            bin_size=args.bin_size,
            bins_per_job=args.j,
            head=args.head)

        for i, result in enumerate(workers.imap_unordered(count_fragments_binned,
                                                          job_gen)):
            counts.update(result)
            print(f'\r{i}/{job_total}', end='')

        if args.o.endswith('.pickle.gz'):
            pd.DataFrame(counts).to_pickle(args.o)
        elif args.o.endswith('.csv'):
            pd.DataFrame(counts).to_csv(args.o)
        else:
            with gzip.open(args.o, 'wb') as out:
                pickle.dump(counts, out)
