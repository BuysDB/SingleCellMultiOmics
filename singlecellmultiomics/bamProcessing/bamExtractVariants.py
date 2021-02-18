#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import singlecellmultiomics
from collections import Counter
from singlecellmultiomics.bamProcessing import sorted_bam_file, has_variant_reads
from singlecellmultiomics.molecule import NlaIIIMolecule,MoleculeIterator,train_consensus_model,get_consensus_training_data, Molecule
from singlecellmultiomics.fragment import NlaIIIFragment, Fragment
import pysam
import collections
import numpy as np
import pandas as pd
import itertools
from singlecellmultiomics.alleleTools import AlleleResolver
import os
from glob import glob
import argparse
from collections import Counter
import multiprocessing
import pickle
import gzip
from contextlib import ExitStack

class VariantWrapper:
    def __init__(self, variant, pos=None,contig=None,ref=None,alts=None,qual=0):
        if pos is not None:
            self.pos = pos
            self.contig = contig
            self.chrom= contig
            self.ref = ref
            self.alts = alts
            self.qual = qual

        else:
            self.pos = variant.pos
            self.contig = variant.contig
            self.chrom= variant.contig
            self.ref = variant.ref
            self.alts = variant.alts
            self.qual = variant.qual
    def __repr__(self):
        return f'VARIANT: {self.contig}:{self.pos} {self.ref} {self.alts} {self.qual}'



def job_gen( induced_variants_path, germline_variants_path,
        germline_variants_sample, alignments_path, block_size = 100, n=None,
        contig=None, completed=None,min_qual=None,germline_bam_path=None,
        MAX_REF_MOLECULES=1000,window_radius=600,max_buffer_size=100_000, debug_bam_folder=None ):
    """
    Job generator

    block_size(int) : variants per block
    n(int) : amount of blocks to generate
    min_qual(float) : minimum quality score of variants to process
    contig: contig to generate jobs for
    completed(set): set of locations which should be skipped
    """

    i=0
    with pysam.VariantFile(induced_variants_path,ignore_truncation=True) as sc_calls:

        vlist = []
        for record in sc_calls:
            if contig is not None and record.chrom!=contig:
                continue

            if completed is not None and  (record.chrom, record.pos) in completed:
                continue

            if min_qual is not None and record.qual<min_qual:
                continue

            if len(record.alts[0])!=1 or len(record.ref)!=1:
                continue

            k = (record.chrom, record.pos)

            vlist.append(VariantWrapper(record))

            if len(vlist)>=block_size:
                #f'./{extraction_folder}/variants_extracted_0_NLA_{i}.bam'
                yield (vlist, alignments_path, None, 'NLA', germline_variants_path,
                    germline_variants_sample, germline_bam_path,
                    window_radius, MAX_REF_MOLECULES,max_buffer_size, debug_bam_folder)

                vlist = []
                i+=1
                if n is not None and i>=n:
                    break
        if len(vlist):
            yield (vlist, alignments_path, None, 'NLA', germline_variants_path,
                germline_variants_sample, germline_bam_path,
                window_radius, MAX_REF_MOLECULES,max_buffer_size, debug_bam_folder)


def get_molecule_base_calls(molecule, variant):
    c = molecule.get_consensus()

    if not (variant.chrom, variant.pos-1) in c:
        return None
    if c[(variant.chrom, variant.pos-1)]==variant.ref:
        return variant.ref, molecule.get_mean_base_quality(variant.chrom, variant.pos-1, variant.ref)
    elif c[(variant.chrom, variant.pos-1)]==variant.alts[0]:
        return variant.alts[0], molecule.get_mean_base_quality(variant.chrom, variant.pos-1, variant.alts[0])

def get_phased_variants(molecule,resolver=None):

    if resolver is None:
        resolver = molecule.allele_resolver

    haplotype = molecule.get_allele(
            return_allele_informative_base_dict=True,
            allele_resolver=resolver)

    return [(chromosome, position, base)
            for allele, bps in haplotype.items()
            for chromosome, position, base in bps]


def filter_alt_calls(alt_phased: collections.Counter, threshold: float):
    """
    Filter the counter alt-phased
    """
    total_per_pos = Counter()
    for (phasedchrom, phased_pos, phased_base),obs in alt_phased.most_common():
        total_per_pos[(phasedchrom, phased_pos)] += obs

    return [(phasedchrom, phased_pos, phased_base)
    for (phasedchrom, phased_pos, phased_base),obs in alt_phased.most_common()
     if obs/total_per_pos[(phasedchrom, phased_pos)] >= threshold
    ]



def recall_variants(args):

    variants, alignment_file_path, target_path, mode, germline_variants_path, germline_variants_sample, germline_bam_path, window_radius, MAX_REF_MOLECULES,max_buffer_size, debug_bam_folder = args

    window_radius = 600
    MAX_REF_MOLECULES = 5_000  # Maximum amount of reference molecules to process.
    # This is capped for regions to which many reads map (mapping artefact)

    variant_calls = dict() # cell->(chrom,pos) +/- ?
    phased_variants = dict()

    ### Set up molecule iterator (1/2)
    if mode== 'NLA':
        mc = NlaIIIMolecule
        fc = NlaIIIFragment
    else:
        mc = Molecule
        fc = Fragment

    ###
    locations_done=set()
    alignments = pysam.AlignmentFile(alignment_file_path,threads=4)
    if germline_bam_path is not None:
        germline_alignments =  pysam.AlignmentFile(germline_bam_path,threads=4)

    for variant in variants:

        # Check if the variant is present in the germline bam file (if supplied)
        if germline_bam_path is not None and has_variant_reads(
                                                        germline_alignments,
                                                        variant.chrom,
                                                        variant.pos-1,
                                                        variant.alts[0],
                                                        min_reads=1,
                                                        stepper='nofilter'):
                #print(f'FOUND IN GERMLINE {variant}')
                continue

        #print(variant)
        overlap = False
        reference_start = max(0, variant.pos - window_radius)
        reference_end = variant.pos + window_radius
        contig = variant.contig

        variant_key = (contig, variant.pos, variant.ref, variant.alts[0] )

        #print(contig,reference_start,reference_end,variant.alts[0],variant.ref)
        ### Set up allele resolver
        unphased_allele_resolver = singlecellmultiomics.alleleTools.AlleleResolver(
            use_cache=False,
            phased=False,
            verbose = True)

        if germline_variants_path is not None:
            with pysam.VariantFile(germline_variants_path) as germline:
                for i, ar_variant in enumerate(germline.fetch(
                        variant.chrom, reference_start, reference_end )):

                    if germline_variants_sample is None:
                        # If any of the samples is not heterozygous: continue
                        if any( (ar_variant.samples[sample].alleles!=2 for sample in ar_variant.samples) ):
                            continue
                    elif len(set(ar_variant.samples[germline_variants_sample].alleles))!=2:
                        continue
                    unphased_allele_resolver.locationToAllele[ar_variant.chrom][ar_variant.pos - 1] = {
                                ar_variant.alleles[0]: {'U'}, ar_variant.alleles[1]: {'V'}
                                }
        ####

        ref_phased = Counter()
        alt_phased = Counter()


        ###

        with ExitStack() as e_stack:

            if debug_bam_folder is not None:
                output_bam = e_stack.enter_context( singlecellmultiomics.bamProcessing.sorted_bam_file(
                    f'{debug_bam_folder}/{"_".join((str(x) for x in variant_key))}.bam', origin_bam=alignments))
            else:
                output_bam = None

            ### Set up molecule iterator (2/2)
            try:
                molecule_iter = MoleculeIterator(
                    alignments,
                    mc,
                    fc,
                    contig=contig,
                    start=reference_start,
                    end=reference_end,
                    molecule_class_args={
                       'allele_resolver':unphased_allele_resolver,
                        'max_associated_fragments':40,
                    },
                    max_buffer_size=max_buffer_size
                )

                reference_called_molecules = [] # molecule, phase

                extracted_base_call_count = 0
                alt_call_count = 0
                for mi,molecule in enumerate(molecule_iter):
                    base_call = get_molecule_base_calls(molecule, variant)
                    if base_call is None:
                        continue
                    extracted_base_call_count+=1
                    base, quality = base_call
                    call = None
                    if base==variant.alts[0]:
                        call='A'
                        alt_call_count+=1
                        if molecule.sample not in variant_calls:
                            variant_calls[molecule.sample] = {}
                        variant_calls[molecule.sample][variant_key] = 1

                    elif base==variant.ref:
                        call='R'

                    if debug_bam_folder is not None:
                        # Write allele-call
                        if call is None:
                            molecule.set_meta('ac','UNK')
                        else:
                            molecule.set_meta('ac', call if call != 'R' else 'UR') # We dont know yet if this is truly
                        # # reference at the allele position or just the uninformative allele

                    if call is None:
                        if output_bam is not None:
                            molecule.write_pysam(output_bam)
                        continue

                    # Obtain all germline variants which are phased :
                    phased = get_phased_variants(molecule, unphased_allele_resolver)

                    if call == 'R' and len(phased) > 0:
                        # If we can phase the alternative allele to a germline variant
                        # the reference calls can indicate absence
                        if len(reference_called_molecules) < MAX_REF_MOLECULES:
                            reference_called_molecules.append((molecule, phased))
                        else:
                            if output_bam is not None:
                                molecule.write_pysam(output_bam)
                    else:
                        if output_bam is not None:
                            molecule.write_pysam(output_bam)

                    for chrom, pos, base in phased:
                        if call == 'A':
                            alt_phased[(chrom, pos, base)] += 1

                        elif call == 'R':
                            ref_phased[(chrom, pos, base)] += 1

            except MemoryError:
                print(f"Buffer exceeded for {variant.contig} {variant.pos}")
                continue

            #print(mi,extracted_base_call_count,alt_call_count)
            if len(alt_phased) > 0 and len(reference_called_molecules):
                # Clean the alt_phased variants for variants which are not >90% the same
                alt_phased_filtered = filter_alt_calls(alt_phased, 0.9)
                #print(alt_phased_filtered)
                phased_variants[variant_key] = alt_phased_filtered
                for molecule, phased_gsnvs in reference_called_molecules:
                    for p in phased_gsnvs:
                        if p in alt_phased_filtered:
                            if not molecule.sample in variant_calls:
                                variant_calls[molecule.sample] = {}
                            variant_calls[molecule.sample][variant_key] = 0

                            if debug_bam_folder is not None:
                                molecule.set_meta('S0', True)
                                molecule.set_meta('ac', 'R')

                            break
                # And write:
                if output_bam is not None:
                    for molecule, phased_gsnvs in reference_called_molecules:
                        molecule.write_pysam(output_bam)


            locations_done.add(variant_key)


    alignments.close()
    return variant_calls, locations_done, phased_variants


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract variants from single cells.
    """)
    argparser.add_argument('bamfile', metavar='bamfile', type=str)

    argparser.add_argument('-extract', help="vcf file with variants to extract", required=True)
    argparser.add_argument('-germline', help="vcf file with germline variants to potentially phase with", required=False)
    argparser.add_argument('-germline_sample', help="germline sample in supplied vcf file")
    argparser.add_argument('-germline_bam', help="germline bam file (no variant reads are allowed in this file)", default=None)
    argparser.add_argument('-debug_bam_folder', help="Path to folder to write informative alignments to")


    argparser.add_argument(
        '-o',
        type=str,
        required=True,
        help='output path, ends in .pickle.gz, or .csv')

    argparser.add_argument(
        '-po',
        type=str,
        required=True,
        help='Phased gsnv output path, ends in .pickle.gz')


    argparser.add_argument('-head', type=int, help='Process only the first N*job_size variants')
    argparser.add_argument('-t', type=int,default=8,help='Threads')
    argparser.add_argument('-minqual', type=float,help='Min variant quality to extract (from the -extract vcf file)')
    argparser.add_argument('-jobsize', type=int,default=5,help='Amount of variants being processed per Thread ')

    args = argparser.parse_args()

    if args.debug_bam_folder is not None and not os.path.exists(args.debug_bam_folder):
        os.makedirs(args.debug_bam_folder)

    assert args.o.endswith('.pickle.gz') or args.o.endswith('.csv')

    variant_calls = collections.defaultdict(dict)
    phased = dict() # variant->((chrom,pos,base),())

    print(f'Initialising {args.t} workers')

    jobs = job_gen( induced_variants_path=args.extract,
            germline_variants_path=args.germline,
            germline_variants_sample=args.germline_sample,
            germline_bam_path=args.germline_bam,
            alignments_path=args.bamfile,
            n=args.head,
            block_size=args.jobsize,
            min_qual=args.minqual,
            debug_bam_folder=args.debug_bam_folder
            )

    if args.t==1:

        def dummy_imap(func, args):
            for arg in args:
                yield func(arg)

        for i,(vc,done, alt_phased) in enumerate(dummy_imap(recall_variants, jobs )):

            for cell, calls in vc.items():
                variant_calls[cell].update(calls)

    else:
        with multiprocessing.Pool( args.t  ) as workers:

            print('Collecting variant calls')
            for i,(vc,done, alt_phased) in enumerate(
                workers.imap_unordered(recall_variants,jobs)):

                for cell, calls in vc.items():
                    variant_calls[cell].update(calls)

                # Write phased dict:
                for key,value in alt_phased.items():
                    print(key, value)
                    phased[key] = value

            print(i)
            if i%25==0:
                print('writing intermediate result')
                df = pd.DataFrame(variant_calls).T.sort_index()
                if args.o.endswith('.csv'):
                    df.to_csv(args.o)
                else:
                    df.to_pickle(args.o)

                with gzip.open(args.po,'wb') as gf:
                    pickle.dump(phased, gf)

    print('Finished collecting variant calls')
    # Write variants to output pickle file:
    print('Writing to output file')
    df = pd.DataFrame(variant_calls).T.sort_index()
    if args.o.endswith('.csv'):
        df.to_csv(args.o)
    else:
        df.to_pickle(args.o)

    with gzip.open(args.po, 'wb') as gf:
        pickle.dump(phased, gf)
    print('Finished')
