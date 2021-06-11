#!/usr/bin/env python
# -*- coding: utf-8 -*-

from singlecellmultiomics.molecule import TAPS, TAPSCHICMolecule, MoleculeIterator
from singlecellmultiomics.fragment import CHICFragment
import singlecellmultiomics.fragment
import singlecellmultiomics.molecule
from singlecellmultiomics.alleleTools import AlleleResolver
from singlecellmultiomics.bamProcessing import get_reference_path_from_bam
import pysam
from glob import glob
import numpy as np
import pandas as pd
from multiprocessing import Pool
from singlecellmultiomics.utils import pool_wrapper
from singlecellmultiomics.bamProcessing.bamBinCounts import blacklisted_binning_contigs
from itertools import product
from pysamiterators import CachedFasta
import pyBigWig
import argparse


def get_methylation_calls(bam, contig, fetch_start, fetch_end, start, stop, mq, bin_size, sample_mapping_function):


    molecule_class_args['reference'] = (
        pysam.FastaFile(get_reference_path_from_bam(bam))
    )

    ar = molecule_class_args.get('allele_resolver')

    with pysam.AlignmentFile(bam) as aln:

        l = stop-start
        nbins = int(l/bin_size)+1
        methylated = dict()
        support = dict()

        for molecule in MoleculeIterator(
                aln,
                molecule_class,
                fragment_class,
                molecule_class_args=molecule_class_args,
                fragment_class_args=fragment_class_args,
            contig=contig,
            start=fetch_start,
            stop=fetch_end

        ):
            if molecule.get_mean_mapping_qual()<mq:
                continue

            sample = sample_mapping_function(molecule.sample)
            if sample is None:
                continue

            if type(sample) is not list and type(sample) is not tuple:
                sample_list = (sample, )
            else:
                sample_list = sample

            for sample in sample_list:
                if not sample in support:
                    support[sample] = np.zeros( nbins )
                    methylated[sample] = np.zeros( nbins )

                for (contig,pos),meta in molecule.methylation_call_dict.items():
                    if meta['context'] not in 'Zz':
                        continue

                    if pos<start or pos>=stop:
                        continue

                    bx = int( (pos-start)/bin_size)
                    support[sample][bx] += 1
                    if meta['context']=='Z':
                        methylated[sample][bx] += 1


    molecule_class_args['reference'].close()
    return contig, start, stop, methylated, support


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract TAPS methylation calls from BAM file')

    argparser.add_argument('alignmentfiles', type=str,nargs='+')
    argparser.add_argument(
        '-o',
        type=str,
        help='prefix of output bigiwig')

    argparser.add_argument(
        '-ref',
        type=str,
        help='path to reference fasta file, auto detected from bamfile')
    argparser.add_argument(
        '--every_fragment_as_molecule',
        action='store_true',
        help='Assign every fragment as a molecule, this effectively disables UMI deduplication and consensus calls based on multiple fragments')
    argparser.add_argument(
        '-min_support',
        type=int,
        help='Minimum amount of molcules for a bin to be taken into account',default=5)

    pseudobulk_gr = argparser.add_argument_group('Pseudobulk settings')
    pseudobulk_gr.add_argument(
        '-pseudobulk_BI_csv',
        type=str,
        help="""Path to a CSV file which contains for every barcode index (BI tag) to what group it belongs.
         The CSV file has no header and two columns, the first column contains the barcode index,
        the second the sample name. Multiple barcode indices can share the same sample name, this will create a pseudobulk signal"""
        )

    pseudobulk_gr.add_argument(
        '-pseudobulk_SM_csv',
        type=str,
        help="""Path to a CSV file which contains for every barcode index (SM tag) to what group it belongs.
         The CSV file has no header and two columns, the first column contains the sample name,
        the second the target sample name. Multiple barcode indices can share the same sample name, this will create a pseudobulk signal"""
        )

    pseudobulk_gr.add_argument(
        '--single_sample',
        action='store_true',
        help="""Each sample/cell in the supplied bam files will get an output bigwig file"""
        )


    bw_gr = argparser.add_argument_group('Bigwig output settings')


    bw_gr.add_argument(
        '-smooth_window',
        type=int,
        help='Smoothing window for smoothed bw',default=10)


    bw_gr.add_argument(
        '-bin_size',
        type=int,
        help='Methylation bin size',default=400)



    argparser.add_argument(
        '-worker_bin_size',
        type=int,
        default=4_000_000,
        help='Size of genomic regions processed per worker')

    argparser.add_argument(
        '-t',
        type=int,
        default=4,
        help='Amount of processes used')

    argparser.add_argument(
        '-features',
        type=str,
        help='Annotate cut locations with these features, path to gtf file')


    argparser.add_argument(
        '-min_phred_score',
        type=int,
        help='Do not call methylation for bases with a phred score lower than "min_phred_score"',default=None)


    argparser.add_argument(
        '-dove_R1_distance',
        type=int,
        help='Do not call methylation N bases from the end of R1',default=8)

    argparser.add_argument(
        '-dove_R2_distance',
        type=int,
        help='Do not call methylation N bases from the end of R2',default=8)

    argparser.add_argument(
        '-skip_last_n_cycles_R1',
        type=int,
        help='Do not call methylation N bases from the end of R1',default=5)



    argparser.add_argument(
        '-skip_first_n_cycles_R1',
        type=int,
        help='Do not call methylation N bases from the start of R1',default=5)

    argparser.add_argument(
        '-skip_last_n_cycles_R2',
        type=int,
        help='Do not call methylation N bases from the end of R2',default=5)

    argparser.add_argument(
        '-skip_first_n_cycles_R2',
        type=int,
        help='Do not call methylation N bases from the start of R2',default=5)


    argparser.add_argument('-minmq', type=int, default=50)
    argparser.add_argument(
        '-contig',
        type=str,
        help='contig to run on, all when not specified')
    argparser.add_argument('-method', type=str, help='nla, or chic')


    argparser.add_argument(
        '--allow_single_end',
        action='store_true',
        help='Allow single end reads')

    allele_gr = argparser.add_argument_group('alleles')
    allele_gr.add_argument('-alleles', type=str, help="Phased allele file (VCF)")
    allele_gr.add_argument(
        '-allele_samples',
        type=str,
        help="Comma separated samples to extract from the VCF file. For example B6,SPRET")
    allele_gr.add_argument(
        '-unphased_alleles',
        type=str,
        help="Unphased allele file (VCF)")
    allele_gr.add_argument(
        '--haplo_molecule_assignment',
        action='store_true',
        help='Take allele information into account during molecule assignment ')

    allele_gr.add_argument(
    '--set_allele_resolver_verbose',
    action='store_true',
    help='Makes the allele resolver print more')
    allele_gr.add_argument(
        '--use_allele_cache',
        action='store_true',
        help='''Write and use a cache file for the allele information. NOTE: THIS IS NOT THREAD SAFE! Meaning you should not use this function on multiple libraries at the same time when the cache files are not yet available.
            Once they are available there is not thread safety issue anymore''')

    args = argparser.parse_args()
    alignments = pysam.AlignmentFile(args.alignmentfiles[0])


    taps=TAPS()

    bams = args.alignmentfiles
    if args.ref is None:
        args.ref = get_reference_path_from_bam(alignments)
        if args.ref is None:
            raise ValueError("Supply reference, -ref")


    fragment_class_args={'umi_hamming_distance': 0,
                         'no_umi_cigar_processing':False}

    molecule_class_args = {

        'taps': taps,
        'taps_strand':'R',
        'min_max_mapping_quality': args.minmq
    }

    if args.method == 'nla':
        molecule_class = singlecellmultiomics.molecule.TAPSNlaIIIMolecule
        fragment_class = singlecellmultiomics.fragment.NlaIIIFragment
        molecule_class_args.update({'site_has_to_be_mapped': True})
    elif args.method == 'chic':
        molecule_class = singlecellmultiomics.molecule.TAPSCHICMolecule
        fragment_class = singlecellmultiomics.fragment.CHICFragment
    else:
        raise ValueError("Supply 'nla' or 'chic' for -method")


    molecule_class_args['methylation_consensus_kwargs'] = {

        'skip_first_n_cycles_R1':args.skip_first_n_cycles_R1,
        'skip_first_n_cycles_R2':args.skip_first_n_cycles_R2,
        'dove_R1_distance':args.dove_R1_distance,
        'dove_R2_distance':args.dove_R2_distance,
        'skip_last_n_cycles_R1':args.skip_last_n_cycles_R1,
        'skip_last_n_cycles_R2':args.skip_last_n_cycles_R2,
        'min_phred_score':args.min_phred_score
        }

    if args.allow_single_end:
        # Single end base calls are "unsafe", allow them :
        molecule_class_args['allow_unsafe_base_calls'] = True
        fragment_class_args['single_end'] = True


    ignore_conversions = set([('C', 'T'), ('G', 'A')])
    if args.alleles is not None and args.alleles!='none':
        molecule_class_args['allele_resolver'] = AlleleResolver(
            args.alleles,
            select_samples=args.allele_samples.split(',') if args.allele_samples is not None else None,
            lazyLoad=True,
            use_cache=args.use_allele_cache,
            verbose = args.set_allele_resolver_verbose,
            ignore_conversions=ignore_conversions)
    ar = molecule_class_args.get('allele_resolver')

    #### start processing:

    bin_size = args.bin_size
    extraction_bin_size= args.worker_bin_size
    met = {}
    sup = {}
    samples = set()

    alias = args.o
    bams = args.alignmentfiles
    contig_whitelist = None
    if args.contig is not None:
        contig_whitelist=[args.contig]


    if args.single_sample:
        def sample_mapping_function(s):
            return s

    elif args.pseudobulk_BI_csv is not None:

        bi_sample_map = {str(bi):str(sample)
            for bi, sample in pd.read_csv(args.pseudobulk_BI_csv,header=None,index_col=0).iloc[:,0].to_dict().items() }
        def sample_mapping_function(s):
            bi = s.split('_')[-1]
            return bi_sample_map.get(bi)
    elif args.pseudobulk_SM_csv is not None:


        sm_sample_map = dict()
        with open(args.pseudobulk_SM_csv) as ip:
            for line in ip:
                parts = line.strip().split()
                sample = parts[0]
                labels = parts[1:]
                sm_sample_map[sample] = labels

        def sample_mapping_function(s):
            return sm_sample_map.get(s)

        """
        sm_sample_map = {str(sm):str(sample)
            for sm, sample in pd.read_csv(args.pseudobulk_SM_csv,header=None,index_col=0).iloc[:,0].to_dict().items() }
        def sample_mapping_function(s):
            return sm_sample_map.get(s)
        """

    else:
        def sample_mapping_function(s):
            return 'bulk'


    with pysam.AlignmentFile(bams[0]) as aln:
        sizes = dict(zip(aln.references,aln.lengths))

    with Pool(args.t) as workers:
        for _contig, _start, _stop,  mr, sr in workers.imap_unordered(
            pool_wrapper,

            [(
                get_methylation_calls,
                {
                    'bam':bam,
                    'start':start,
                    'stop':end,
                    'bin_size':bin_size,
                    'mq':50,
                    'contig':contig,
                    'fetch_start':fetch_start,
                    'fetch_end':fetch_end,
                    'sample_mapping_function':sample_mapping_function,

                }

            ) for bam,(contig,start,end,fetch_start,fetch_end) in product( bams,
                                                                   list(
                                                                       blacklisted_binning_contigs(
                                                                           bams[0],
                                                                           bin_size=extraction_bin_size,
                                                                           fragment_size=800,
                                                                           contig_whitelist=contig_whitelist
                                                                       )))]

        ):

            start_index = int(_start/bin_size)
            l = _stop-_start
            nbins = int(l/bin_size)+1
            stop_index=start_index+nbins


            if not _contig in met:
                met[_contig] = {} # contig-> cell ->numpy array
                sup[_contig] = {}



            for cell, cmr in mr.items():
                if not cell in met[_contig]:
                    met[_contig][cell] = np.zeros( int(sizes[_contig]/bin_size)+1 )
                met[_contig][cell][start_index:stop_index] += cmr


            for cell, smr in sr.items():
                if not cell in sup[_contig]:
                    sup[_contig][cell] = np.zeros( int(sizes[_contig]/bin_size)+1 )
                sup[_contig][cell][start_index:stop_index] += smr
                samples.add(cell)
    ### write to output:
    # Write to the bigwig file:
    min_support=  args.min_support
    smooth_window = args.smooth_window

    for sample in samples:
        with pyBigWig.open(f'{alias}_{sample}.bw','w') as out, \
             pyBigWig.open(f'{alias}_support_{sample}.bw','w') as supp_out, \
             pyBigWig.open(f'{alias}_smooth_{sample}.bw','w') as outs:

            out.addHeader(list(zip(sizes.keys(), sizes.values())))
            outs.addHeader(list(zip(sizes.keys(), sizes.values())))
            supp_out.addHeader(list(zip(sizes.keys(), sizes.values())))

            for contig in sizes:
                if not contig in sup:
                    continue
                if not sample in sup[contig]:
                    continue

                with np.errstate(invalid='ignore'):
                    y = met[contig][sample] / sup[contig][sample]
                coordinates = np.arange(len(y)*bin_size,step=bin_size)


                sel = sup[contig][sample]>min_support
                coordinates_defined = coordinates[sel]
                y=  y[sel]

                if len(y)<=1:
                    continue

                out.addEntries(
                            [contig]*len(y), #Contig
                            list(coordinates_defined), #Start
                            ends= list(coordinates_defined+bin_size), #end
                            values= list(y))

                # Write support:
                supp_out.addEntries(
                            [contig]*len(y), #Contig
                            list(coordinates_defined), #Start
                            ends= list(coordinates_defined+bin_size), #end
                            values= list(sup[contig][sample][sel]))


                ser = pd.Series( np.interp(coordinates, xp=coordinates_defined, fp=y) ,
                     index=coordinates
                     ).rolling(window=smooth_window,center=True).mean()
                ser = ser[ser>=0]
                if len(ser)>1:
                    outs.addEntries(
                            [contig]*len(ser), #Contig
                            list(ser.index), #Start
                            ends= list(ser.index+bin_size), #end
                            values= list(ser.values))
