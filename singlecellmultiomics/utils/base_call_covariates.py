#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import Pool
from pysam import AlignmentFile, FastaFile
import pysam
from singlecellmultiomics.bamProcessing.bamBinCounts import blacklisted_binning_contigs
from singlecellmultiomics.utils.sequtils import reverse_complement, get_context
from singlecellmultiomics.utils import prob_to_phred
from pysamiterators import CachedFasta
from array import array
from uuid import uuid4
from singlecellmultiomics.bamProcessing import merge_bams, get_contigs_with_reads, has_variant_reads
import argparse
import pickle
import gzip
import pandas as pd
import numpy as np
import os

def get_covariate_key(read, qpos, refpos, reference, refbase, cycle_bin_size=3, k_rad=1):

    qual = read.query_qualities[qpos]
    qbase = read.query_sequence[qpos]

    if qbase == 'N':
        return None

    context = get_context(read.reference_name, refpos, reference, qbase, k_rad)
    if 'N' in context or len(context)!=(2*k_rad+1):
        context=get_context(read.reference_name, refpos, reference, qbase, 0)
        assert len(context)==1

    if 'N' in context:
        return None


    if read.is_reverse:
        cycle = len(read.query_sequence) - qpos  # -1
        context = reverse_complement(context)
    else:
        cycle = qpos + 1

    return qual, read.is_read2, int(round(cycle / cycle_bin_size)) * cycle_bin_size, context


def add_keys_excluding_context_covar(covariates, covar_phreds, k_rad=1):
    k_rad = 1
    for base in 'ACTG':
        pd.DataFrame({ (key[0],key[1],key[2],key[3][k_rad]):value for key,value in covariates.items() if len(key[3])==k_rad*2+1 and key[3][k_rad]==base }).T

    index = []
    values = []
    for key,(t,f) in [((key[0],key[1],key[2],key[3][k_rad]), value) for key,value in covariates.items() if len(key[3])==k_rad*2+1 and key[3][k_rad]==base ]:
        index.append(key)
        values.append([t,f])

    for indices, base_data in pd.DataFrame(values,index=pd.MultiIndex.from_tuples(index)).groupby(level=(0,1,2,3)):
        t,f = base_data.sum()
        covar_phreds[indices] = prob_to_phred( f/(t+f) )

def covariate_obs_to_phreds(covariates, k_rad):
    covar_phreds = {}

    for k,(p_base_true,p_base_false) in covariates.items():
        covar_phreds[k] = prob_to_phred( p_base_false/(p_base_true+p_base_false) )

    # Add values for context not available:
    add_keys_excluding_context_covar(covariates, covar_phreds,k_rad=k_rad)

    # Add value for when everything fails (just mean):
    t,f = pd.DataFrame(covariates).mean(1)
    covar_phreds[None] = prob_to_phred( f/(t+f) )

    return covar_phreds

def extract_covariates(bam_path: str,
                       reference_path: str,
                       contig: str,
                       start: int,
                       end: int,
                       start_fetch: int,
                       end_fetch: int,
                       filter_kwargs: dict,
                       covariate_kwargs: dict):
    """
    Count mismatches and matches for similar base-calls

    Returns:
        match_mismatch(dict) : dictionary ( covariate_key: [mismatches, matches], .. )
    """
    # known is a set() containing locations of known variation (snps)
    # @todo: extend to indels
    global known  # <- Locations, set of (contig, position) tuples to ignore

    joined = dict()

    # Filters which select which reads are used to estimate covariates:
    min_mapping_quality = filter_kwargs.get('min_mapping_quality', 0)
    deduplicate = filter_kwargs.get('deduplicate', False)
    filter_qcfailed = filter_kwargs.get('filter_qcfailed', False)
    variant_blacklist_vcf_files = filter_kwargs.get('variant_blacklist_vcf_files', None)

    # Obtain all variants in the selected range:
    blacklist = set()
    if variant_blacklist_vcf_files is not None:

        for path in variant_blacklist_vcf_files:
            with pysam.VariantFile(path) as bf:
                for record in bf.fetch(contig, start_fetch, end_fetch):
                    blacklist.add(record.pos)


    with AlignmentFile(bam_path) as alignments,  FastaFile(reference_path) as fa:
        reference = CachedFasta(fa)  # @todo: prefetch selected region
        for read in alignments.fetch(contig, start_fetch, end_fetch):
            if (deduplicate and read.is_duplicate) or \
                    (read.is_qcfail and filter_qcfailed) or \
                    (read.mapping_quality < min_mapping_quality):
                continue

            for qpos, refpos, refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):

                if refpos > end or refpos < start:  # Prevent the same location to be counted multiple times
                    continue

                if refpos in blacklist:
                    continue

                refbase = refbase.upper()
                if refbase == 'N' or (read.reference_name, refpos) in known:
                    continue

                key = get_covariate_key(read, qpos, refpos, reference, refbase, **covariate_kwargs)
                if key is None:
                    continue

                matched = (refbase == read.query_sequence[qpos])
                try:
                    joined[key][matched] += 1
                except KeyError:
                    if matched:
                        joined[key] = array('l', [0, 1])
                    else:
                        joined[key] = array('l', [1, 0])
    return joined


def extract_covariates_wrapper(kwargs):
    return extract_covariates(**kwargs)


def extract_covariates_from_bam(bam_path, reference_path, known_variants, n_processes=None, bin_size=10_000_000,
            min_mapping_quality = 40,
            deduplicate = True,
            filter_qcfailed = True,
            variant_blacklist_vcf_files = None
            ):

    global known
    known = known_variants

    joined = dict()

    job_generation_args = {
        'contig_length_resource': bam_path,
        'bin_size': bin_size,
        'fragment_size': 0}

    filter_kwargs = {
        'min_mapping_quality': 40,
        'deduplicate': True,
        'filter_qcfailed': True,
        'variant_blacklist_vcf_files':variant_blacklist_vcf_files
    }

    covariate_kwargs = {
        'cycle_bin_size': 3,
        'k_rad' : 1
    }

    jobs_total = sum(1 for _ in (blacklisted_binning_contigs(**job_generation_args)))

    with Pool(n_processes) as workers:

        for i, r in enumerate(
                        workers.imap_unordered(extract_covariates_wrapper, (
                                   {
                                         'bam_path': bam_path,
                                         'reference_path': reference_path,
                                         'contig': contig,
                                         'start': start,
                                         'end': end,
                                         'start_fetch': start_fetch,
                                         'end_fetch': end_fetch,
                                         'filter_kwargs': filter_kwargs,
                                         'covariate_kwargs': covariate_kwargs
                                    }
                                   for contig, start, end, start_fetch, end_fetch in
                                             blacklisted_binning_contigs(**job_generation_args)))):
            print(round(100 * (i / jobs_total), 1), end='\r')

            for key, tf in r.items():
                try:
                    joined[key][0] += tf[0]
                    joined[key][1] += tf[1]
                except KeyError:
                    joined[key] = array('l', [0, 0])
                    joined[key][0] += tf[0]
                    joined[key][1] += tf[1]
    return joined


def recalibrate_base_calls(read, reference, joined_prob, covariate_kwargs):
    # @todo: make copy to save phred scores of soft-clipped bases
    # This array will contain all recalibrated phred scores:
    new_qualities = array('B', [0] * len(read.query_qualities))

    # Iterate all aligned pairs and replace phred score:

    for qpos, refpos in read.get_aligned_pairs(matches_only=True, with_seq=False):

        key = get_covariate_key(read, qpos, refpos, reference, None, **covariate_kwargs)
        try:
            phred = joined_prob[key]
        except KeyError:
            phred = joined_prob[None]
        new_qualities[qpos] = phred

    read.query_qualities = new_qualities


def _recalibrate_reads(bam_path, reference_path, contig, start, end, covariate_kwargs, **kwargs):
    # Recalibrate the reads in bam_path

    global joined_prob  # Global to share over multiprocessing
    # joined_prob contains  P(error| d), where d is a descriptor generated by get_covariate_key

    o_path = f'out_{uuid4()}.bam'

    # Open source bam file:
    with AlignmentFile(bam_path) as alignments, FastaFile(reference_path) as fa:
        # @todo: extract only selected region from fasta file:
        reference = CachedFasta(fa)
        # Open target bam file:
        with AlignmentFile(o_path, header=alignments.header, mode='wb') as out:
            # Iterate all reads in the source bam file:
            for read in alignments.fetch(contig, start, end):
                recalibrate_base_calls(read, reference, joined_prob, covariate_kwargs)
                out.write(read)

    pysam.index(o_path)
    return o_path


def __recalibrate_reads(kwargs):
    return _recalibrate_reads(**kwargs)


def recalibrate_reads(bam_path, target_bam_path, reference_path, n_processes, covariates,  covariate_kwargs, intermediate_bam_size=20_000_000):
    job_generation_args = {
        'contig_length_resource': bam_path,
        'bin_size': intermediate_bam_size,
        'fragment_size': 0
    }
    global joined_prob
    joined_prob = covariates

    print(len(covariates), 'discrete elements')
    with Pool(n_processes) as workers:
        intermediate_bams = list( workers.imap_unordered(__recalibrate_reads, (
                        {
                            'bam_path': bam_path,
                            'reference_path': reference_path,
                            'contig': contig,
                            'start': None,
                            'end': None,
                           # 'start_fetch': start_fetch,
                            #'end_fetch': end_fetch,
                            'covariate_kwargs': covariate_kwargs
                        }
                        for contig in list(get_contigs_with_reads(bam_path)))))
                        #for contig, start, end, start_fetch, end_fetch in
                        #blacklisted_binning_contigs(**job_generation_args))))

        merge_bams(intermediate_bams, target_bam_path)




if __name__=='__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Obtain base calling biases using complete posterior distribution and perform corrections. Both identifying the covariates and the recallibration are multithreaded""")
    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('-reference', help="Path to reference fasta file used to generate the bamfile", required=True)
    argparser.add_argument('-known', help="vcf file with known variation", required=False)


    argparser.add_argument('-threads', help="Amount of threads to use. Uses all when not set")

    argparser.add_argument(
        '-covariates_out',
        type=str,
        help='Write covariates to this file, ends in .pickle.gz. ')

    argparser.add_argument(
        '-covariates_in',
        type=str,
        help='Read in existing covariates, ends in .pickle.gz')

    argparser.add_argument(
        '-bam_out',
        type=str,
        required=False,
        help='Write corrected bam file here')

    argparser.add_argument(
        '--f',
        action= 'store_true')

    args = argparser.parse_args()

    assert args.bam_out!=args.bamfile, 'The input bam file name cannot match the output bam file'

    if args.covariates_in is None and args.covariates_out is not None and args.bam_out is None and os.path.exists(args.covariates_out) and not args.f:
        print('Output covariate file already exists. Use -f to overwrite.')
        exit()

    # Set defaults  when nothing is supplied
    if args.covariates_in is None and args.bam_out is None and args.covariates_out is None:
        args.bam_out = args.bamfile.replace('.bam','.recall.bam')
        if args.covariates_out is None:
            args.covariates_out = args.bamfile.replace('.bam','.covariates.pickle.gz')


    if args.covariates_in is not None:
        print(f'Loading covariates from {args.covariates_in} ')
        with gzip.open(args.covariates_in,'rb') as i:
            covariates = pickle.load(i)
    else:
        covariates = extract_covariates_from_bam(args.bamfile,
                            reference_path=args.reference,
                               variant_blacklist_vcf_files= [] if args.known is None else [args.known],
                                         known_variants=set()
                           )
        if args.covariates_out is not None:

            print(f'Writing covariates to {args.covariates_out}')
            with gzip.open(args.covariates_out,'wb') as o:
                pickle.dump(covariates, o)

    if args.bam_out is not None:
        # Create cov probs

        covar_phreds =  covariate_obs_to_phreds(covariates ,k_rad=1)

        print(f'Writing corrected bam file to {args.bam_out}')
        recalibrate_reads(bam_path=args.bamfile,
                            target_bam_path=args.bam_out,
                            reference_path=args.reference,
                            n_processes=args.threads,
                            covariates=covar_phreds,
                            covariate_kwargs = {
                                'cycle_bin_size': 3,
                                'k_rad' : 1
                            })
