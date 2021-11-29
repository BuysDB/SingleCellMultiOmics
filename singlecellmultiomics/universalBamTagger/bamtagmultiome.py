#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
from singlecellmultiomics.molecule import MoleculeIterator, ReadIterator
import singlecellmultiomics
import singlecellmultiomics.molecule
from singlecellmultiomics.molecule.consensus import calculate_consensus
import singlecellmultiomics.fragment
from singlecellmultiomics.bamProcessing.bamFunctions import sorted_bam_file, get_reference_from_pysam_alignmentFile, write_program_tag, MapabilityReader, verify_and_fix_bam,add_blacklisted_region
from singlecellmultiomics.utils import is_main_chromosome
from singlecellmultiomics.utils.submission import submit_job
import singlecellmultiomics.alleleTools
from singlecellmultiomics.universalBamTagger.customreads import CustomAssingmentQueryNameFlagger, BulkFlagger
import singlecellmultiomics.features
from pysamiterators import MatePairIteratorIncludingNonProper, MatePairIterator
from singlecellmultiomics.universalBamTagger.tagging import generate_tasks, prefetch, run_tagging_tasks, write_job_gen_to_bed
from singlecellmultiomics.bamProcessing.bamBinCounts import blacklisted_binning_contigs,get_bins_from_bed_iter
from singlecellmultiomics.utils.binning import bp_chunked
from singlecellmultiomics.bamProcessing import merge_bams, get_contigs_with_reads, sam_to_bam
from singlecellmultiomics.fastaProcessing import CachedFastaNoHandle
from singlecellmultiomics.utils.prefetch import UnitialisedClass
from multiprocessing import Pool
from typing import Generator
import argparse
import uuid
import os
import sys
import colorama
import pkg_resources
import pickle
from datetime import datetime
from time import sleep
# Supress [E::idx_find_and_load] Could not retrieve index file for, see https://github.com/pysam-developers/pysam/issues/939
pysam.set_verbosity(0)

argparser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Assign molecules, set sample tags, set alleles')
argparser.add_argument('bamin', type=str, help='Input BAM file, SAM files can also be supplied but will be converted to BAM')
argparser.add_argument(
    '-ref',
    type=str,
    default=None,
    help="Path to reference fast (autodected if not supplied)")

argparser.add_argument('-o', type=str, help="output bam file", required=True)
argparser.add_argument(
    '-method',
    type=str,
    default=None,
    required=True,
    help="""Protocol to tag, select from:nla, qflag, chic, nla_transcriptome, vasa, cs, cs_feature_counts,  nla_taps ,chic_taps, nla_no_overhang. nla (Data with digested by Nla III enzyme)
    nla (Data with digested by Nla III enzyme)
    qflag (Only add basic tags like sampple and UMI, no molecule assignment)
    chic (Data digested using mnase fusion)
    nla_transcriptome (Data with transcriptome and genome digested by Nla III )
    vasa (VASA transcriptomic data)
    cs (CELseq data, 1 and 2)
    cs_feature_counts (Single end, deduplicate using a bam file tagged using featurecounts, deduplicates a umi per gene)
    fl_feature_counts (deduplicate using a bam file tagged using featurecounts, deduplicates based on fragment position)
    nla_taps (Data with digested by Nla III enzyme and methylation converted by TAPS)
    chic_taps (Data with digested by mnase enzyme and methylation converted by TAPS)
    nla_tapsp_transcriptome (Add feature annotation to nla_ptaps mode )
    nla_taps_transcriptome  (Add feature annotation to nla_taps mode )
    nla_no_overhang (Data with digested by Nla III enzyme, without the CATG present in the reads)
    scartrace (Lineage tracing )
    """)
argparser.add_argument(
    '-qflagger',
    type=str,
    default=None,
    help="Query flagging algorithm, use -bulk to set the same sample and umi to all reads")
argparser.add_argument(
    '--ignore_bam_issues',
    action='store_true',
    help='Ignore truncation')
argparser.add_argument('-custom_flags', type=str, default="MI,RX,bi,SM")

r_select = argparser.add_argument_group('Region selection')
r_select.add_argument('-head', type=int)
r_select.add_argument('-contig', type=str, help='Contig to only process')
r_select.add_argument('-skip_contig', type=str, help='Contigs not to process')
r_select.add_argument('-region_start', type=int, help='Zero based start coordinate of region to process')
r_select.add_argument('-region_end', type=int, help='Zero based end coordinate of region to process')
r_select.add_argument('-blacklist', type=str, help='BED file containing regions to skip')

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
    '--set_allele_resolver_verbose',
    action='store_true',
    help='Makes the allele resolver print more')

allele_gr.add_argument(
    '--haplo_molecule_assignment',
    action='store_true',
    help='Take allele information into account during molecule assignment ')


allele_gr.add_argument(
    '--use_allele_cache',
    action='store_true',
    help='''Write and use a cache file for the allele information. NOTE: THIS IS NOT THREAD SAFE! Meaning you should not use this function on multiple libraries at the same time when the cache files are not yet available.
        Once they are available there is not thread safety issue anymore''')

argparser.add_argument('-molecule_iterator_verbosity_interval',type=int,default=None,help='Molecule iterator information interval in seconds')
argparser.add_argument('--molecule_iterator_verbose', action='store_true', help='Show progress indication on command line')
argparser.add_argument('-stats_file_path',type=str,default=None,help='Path to logging file, ends with ".tsv"')
argparser.add_argument('-jobbed',type=str,default=None,help='Path to location to write multiprocessing region log file')

argparser.add_argument(
    '--multiprocess',
    action='store_true',
    help="Use multiple the CPUs of you system to achieve (much) faster tagging")


argparser.add_argument(
    '--one_contig_per_process',
    action='store_true',
    help="Do not split contigs/chromosomes into chunks for parallel processing. Use this when you want to correctly deduplicate mates which are mapping very far apart. (>50kb)")



argparser.add_argument(
    '-tagthreads',
    type=int,
    help='Amount of processes to use for tagging (--multiprocess needs to be enabled). Uses all available CPUs when not set.'
    )

argparser.add_argument(
    '-max_time_per_segment',
    default=None,
    type=float,
    help='Maximum time spent on a single genomic location')

argparser.add_argument(
    '-temp_folder',
    default='./',
    help="Temp folder location")


fragment_settings = argparser.add_argument_group('Fragment settings')
fragment_settings.add_argument('-read_group_format', type=int, default=0, help="0: Every cell/sequencing unit gets a read group, 1: Every library/sequencing unit gets a read group")
fragment_settings.add_argument('-max_fragment_size', type=int, default=None, help='Reject fragments with a fragment size higher the specified value')
fragment_settings.add_argument(
    '--resolve_unproperly_paired_reads',
    action='store_true',
    help='When enabled bamtagmultiome will look through the complete bam file in a hunt for the mate, the two mates will always end up in 1 molecule if both present in the bam file. This also works when the is_proper_pair bit is not set. Use this option when you want to find the breakpoints of genomic re-arrangements.')
fragment_settings.add_argument(
    '--allow_cycle_shift',
    action='store_true',
    help='NlaIII: Allow fragments to be shifted slightly around the cut site.')

fragment_settings.add_argument(
    '-assignment_radius',
    type=int,
    help='Molecule assignment radius')

fragment_settings.add_argument(
    '-libname',type=str,
    help='Overwrite library name with this value')


fragment_settings.add_argument(
    '--no_rejects',
    action='store_true',
    help='Do not write rejected reads to output file')

fragment_settings.add_argument(
    '--no_overflow',
    action='store_true',
    help='Do not write overflow reads to output file. Overflow reads are reads which are discarded because the molecule reached the maximum capacity of associated fragments')

fragment_settings.add_argument(
    '--no_restriction_motif_check',
    action='store_true',
    help='Do not check for restriction motifs (NLAIII)')


molecule_settings = argparser.add_argument_group('Molecule settings')
molecule_settings.add_argument(
    '-mapfile',
    type=str,
    help='Path to *.safe.bgzf file, used to decide if molecules are uniquely mappable, generate one using createMapabilityIndex.py ')
molecule_settings.add_argument('-umi_hamming_distance', type=int, default=1)
molecule_settings.add_argument(
    '-annotmethod',
    type=int,
    default=1,
    help="Annotation resolving method. 0: molecule consensus aligned blocks. 1: per read per aligned base")

molecule_settings.add_argument(
    '--feature_container_verbose',
        action='store_true',
        help='Make the feature container print more')

cluster = argparser.add_argument_group('cluster execution')
cluster.add_argument(
    '--cluster',
    action='store_true',
    help='split by chromosomes and submit the job on cluster')

cluster.add_argument(
    '-mem',
    default=40,
    type=int,
    help='Memory requested per job')
cluster.add_argument(
    '-time',
    default=52,
    type=int,
    help='Time requested per job')

cluster.add_argument(
    '-sched',
    default='slurm',
    type=str,
    help='Scheduler to use: sge, slurm, local')

cluster.add_argument(
    '-clusterdir',
    type=str,
    default=None,
    help='Folder to store cluster files in (scripts and sterr/out, when not specified a "cluster" folder will be made in same directory as -o')

scartrace_settings = argparser.add_argument_group('Scartrace specific settings')
scartrace_settings.add_argument('-scartrace_r1_primers', type=str, default='CCTTGAACTTCTGGTTGTAG', help='Comma separated list of R1 primers used. Only fragments starting with this read are taken into account.')

tr = argparser.add_argument_group('transcriptome specific settings')
tr.add_argument('-exons', type=str, help='Exon GTF file')
tr.add_argument(
    '-introns',
    type=str,
    help='Intron GTF file, use exonGTFtoIntronGTF.py to create this file')

cg = argparser.add_argument_group('molecule consensus specific settings')
cg.add_argument(
    '--consensus',
    action='store_true',
    help='Calculate molecule consensus read, this feature is _VERY_ experimental')

cg.add_argument('--no_source_reads', action='store_true',
                help='Do not write original reads, only consensus ')



ma = argparser.add_argument_group('Molecule assignment settings')

ma.add_argument(
    '--every_fragment_as_molecule',
    action='store_true',
    help='Assign every fragment as a molecule, this effectively disables UMI deduplication')

ma.add_argument(
    '--no_umi_cigar_processing',
    action='store_true',
    help='Do not use the alignment during deduplication')
ma.add_argument('-max_associated_fragments',type=int, default=None, help="Limit the maximum amount of reads associated to a single molecule.")


def tag_multiome_multi_processing(
        input_bam_path: str,
        out_bam_path: str,
        molecule_iterator: Generator = None, # @todo add custom molecule iterator?
        molecule_iterator_args: dict = None,
        ignore_bam_issues: bool = False,  # @todo add ignore_bam_issues
        head: int = None,  # @todo add head
        no_source_reads: bool = False,  # @todo add no_source_reads
        # One extra parameter is the fragment size:
        fragment_size: int = None,
        # And the blacklist is optional:
        blacklist_path: str = None,
        bp_per_job: int = None,
        bp_per_segment: int = None,
        temp_folder_root: str = '/tmp/scmo',
        max_time_per_segment: int = None,
        use_pool: bool = True,
        one_contig_per_process: bool =False,
        additional_args: dict = None,
        n_threads=None,
        job_bed_file: str = None # Writes job blocks to a bed file for inspection
    ):

    assert bp_per_job is not None
    assert fragment_size is not None
    assert bp_per_segment is not None

    if not os.path.exists(temp_folder_root):
        raise ValueError(f'The path {temp_folder_root} does not exist')
    temp_folder = os.path.join( temp_folder_root , f'scmo_{uuid.uuid4()}' )
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder, exist_ok=True)

    if molecule_iterator_args.get('skip_contigs', None) is not None:
        contig_blacklist = molecule_iterator_args.get('skip_contigs')
    else:
        contig_blacklist = []

    if molecule_iterator_args.get('contig',None) is not None:
        assert molecule_iterator_args.get('start') is None, 'regions are not implemented'
        contig_whitelist = [molecule_iterator_args['contig']]
    else:
        contig_whitelist = [contig for contig in get_contigs_with_reads(input_bam_path) if not contig in contig_blacklist]

    for prune in ['start', 'end', 'contig', 'progress_callback_function','alignments']:
        if prune in molecule_iterator_args:
            del molecule_iterator_args[prune]

    iteration_args = {
        'molecule_iterator_args': molecule_iterator_args,
        'molecule_iterator_class': MoleculeIterator
    }

    # Define the regions to be processed and group into segments to perform tagging on
    if one_contig_per_process:

        #job_gen = [ [(contig,0,contig_len,0,contig_len),] for contig,contig_len in get_contigs_with_reads(input_bam_path, True)  ]
        if blacklist_path is not None:
            raise NotImplementedError('A blacklist is currently incompatible with --multiprocessing in single contig mode')
        job_gen = [ [(contig,None,None,None,None),] for contig,contig_len in get_contigs_with_reads(input_bam_path, True) if contig!='*' ]


    else:
        regions = blacklisted_binning_contigs(
                contig_length_resource=input_bam_path,
                bin_size=bp_per_segment,
                fragment_size=fragment_size,
                blacklist_path=blacklist_path,
                contig_whitelist=contig_whitelist
            )

        # Chunk into jobs of roughly equal size: (A single job will process multiple segments)
        job_gen = bp_chunked(regions, bp_per_job)


    if job_bed_file is not None:

        assert not one_contig_per_process, 'Cannot write bed file with jobs. Each contig is a separate job, there are no bins.'
        job_gen = list(job_gen) # Convert the job generator to a list to allow it to be traversed multiple times
        write_job_gen_to_bed( job_gen, job_bed_file)
        print(f'Wrote job bed file to {job_bed_file}')

    tasks = generate_tasks(input_bam_path=input_bam_path,
                           job_gen=job_gen,
                           iteration_args=iteration_args,
                           temp_folder=temp_folder,
                           additional_args= additional_args,
                           max_time_per_segment=max_time_per_segment)

    # Create header bam:
    temp_header_bam_path = f'{temp_folder}/{uuid.uuid4()}_header.bam'
    with pysam.AlignmentFile(input_bam_path) as input_bam:
        input_header = input_bam.header.as_dict()

        # Write provenance information to BAM header
        write_program_tag(
            input_header,
            program_name='bamtagmultiome',
            command_line=" ".join(
                sys.argv),
            version=singlecellmultiomics.__version__,
            description=f'SingleCellMultiOmics molecule processing, executed at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')


        if blacklist_path is not None:
            for contig, start, end in get_bins_from_bed_iter(blacklist_path):
                add_blacklisted_region(input_header,contig,start,end,source=blacklist_path)

        # Prefetch the genomic resources with the defined genomic interval reducing I/O load during processing of the region
        # this is now done automatically, by tagging.run_tagging_task

        # @todo : Obtain auto blacklisted regions if applicable
        # @todo : Progress indication

        bam_files_generated = []

        if use_pool:
            workers = Pool(n_threads)
            job_generator = workers.imap_unordered(run_tagging_tasks, tasks)

        else:
            job_generator = (run_tagging_tasks(task) for task in tasks)

        total_processed_molecules = 0
        for bam, meta in job_generator:
            if bam is not None:
                bam_files_generated.append(bam)
            if len(meta):
                total_processed_molecules+=meta['total_molecules']
                timeouts = meta.get('timeout_tasks',[])
                for timeout in timeouts:
                    print('blacklisted', timeout['contig'], timeout['start'], timeout['end'])
                    add_blacklisted_region(input_header,
                        contig=timeout['contig'],
                        start=timeout['start'],
                        end=timeout['end']
                    )
            if head is not None and total_processed_molecules>head:
                print('Head was supplied, stopping')
                break
        tagged_bam_generator = [temp_header_bam_path] +  bam_files_generated

        with pysam.AlignmentFile(temp_header_bam_path, 'wb', header=input_header) as out:
            pass

        pysam.index(temp_header_bam_path)
    # merge the results and clean up:
    print('Merging final bam files')
    merge_bams(list(tagged_bam_generator), out_bam_path)
    if use_pool:
        workers.close()

    # Remove the temp dir:
    sleep(5)
    try:
        os.rmdir(temp_folder)
    except Exception:
        sys.stderr.write(f'Failed to remove {temp_folder}\n')

    write_status(out_bam_path, 'Reached end. All ok!')


def tag_multiome_single_thread(
        input_bam_path,
        out_bam_path,
        molecule_iterator = None,
        molecule_iterator_args = None,
        consensus_model = None,
        consensus_model_args={}, # Clearly the consensus model class and arguments should be part of molecule
        ignore_bam_issues=False,
        head=None,
        no_source_reads=False
        ):

    input_bam = pysam.AlignmentFile(input_bam_path, "rb", ignore_truncation=ignore_bam_issues, threads=4)
    input_header = input_bam.header.as_dict()


    # Write provenance information to BAM header
    write_program_tag(
        input_header,
        program_name='bamtagmultiome',
        command_line=" ".join(
            sys.argv),
        version=singlecellmultiomics.__version__,
        description=f'SingleCellMultiOmics molecule processing, executed at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

    print(f'Started writing to {out_bam_path}')

    # de-prefetch all:


    # contig, start, end, start, end , args

    molecule_iterator_args = prefetch(molecule_iterator_args['contig'],
                                    molecule_iterator_args['start'],
                                    molecule_iterator_args['end'],
                                    molecule_iterator_args['start'],
                                    molecule_iterator_args['end'],
                                    molecule_iterator_args)


    molecule_iterator_exec = molecule_iterator(input_bam, **{k:v for k, v in molecule_iterator_args.items()
                                                            if k != 'alignments'})

    print('Params:',molecule_iterator_args)
    read_groups = dict()  # Store unique read groups in this dict



    with sorted_bam_file(out_bam_path, header=input_header, read_groups=read_groups) as out:
        try:
            for i, molecule in enumerate(molecule_iterator_exec):

                # Stop when enough molecules are processed
                if head is not None and (i - 1) >= head:
                    break

                # set unique molecule identifier
                molecule.set_meta('mi', f'{molecule.get_a_reference_id()}_{i}')

                # Write tag values
                molecule.write_tags()

                """
                if unphased_allele_resolver is not None:  # write unphased allele tag:
                    molecule.write_allele_phasing_information_tag(
                        unphased_allele_resolver, 'ua')
                """

                # Update read groups
                for fragment in molecule:
                    rgid = fragment.get_read_group()
                    if not rgid in read_groups:
                        read_groups[rgid] = fragment.get_read_group(True)[1]

                # Calculate molecule consensus
                if consensus_model is not None:
                    calculate_consensus(molecule,
                                        consensus_model,
                                        i,
                                        out,
                                        **consensus_model_args)

                # Write the reads to the output file
                if not no_source_reads:
                    molecule.write_pysam(out)
        except Exception as e:
            write_status(out_bam_path,'FAIL, The file is not complete')
            raise e

        # Reached the end of the generator
        write_status(out_bam_path,'Reached end. All ok!')

def write_status(output_path, message):
    status_path = output_path.replace('.bam','.status.txt')
    with open(status_path,'w') as o:
        o.write(message+'\n')


def run_multiome_tagging_cmd(commandline):
    args = argparser.parse_args(commandline)
    run_multiome_tagging(args)


def run_multiome_tagging(args):
    """
    Run multiome tagging adds molecule information

    Arguments:

        bamin (str) : bam file to process, sam files can also be supplied but will be converted

        o(str) : path to output bam file

        method(str): Protocol to tag, select from:nla, qflag, chic, nla_transcriptome, vasa, cs, nla_taps ,chic_taps, nla_no_overhang, scartrace

        qflagger(str): Query flagging algorithm to use, this algorithm extracts UMI and sample information from your reads. When no query flagging algorithm is specified, the `singlecellmultiomics.universalBamTagger.universalBamTagger.QueryNameFlagger` is used

        method(str) : Method name, what kind of molecules need to be extracted. Select from:
            nla (Data with digested by Nla III enzyme)
            qflag (Only add basic tags like sampple and UMI, no molecule assignment)
            chic (Data digested using mnase fusion)
            nla_transcriptome (Data with transcriptome and genome digested by Nla III )
            vasa (VASA transcriptomic data)
            cs (CELseq data, 1 and 2)
            cs_feature_counts (Single end, deduplicate using a bam file tagged using featurecounts, deduplicates a umi per gene)
            fl_feature_counts (deduplicate using a bam file tagged using featurecounts, deduplicates based on fragment position)
            nla_taps (Data with digested by Nla III enzyme and methylation converted by TAPS)
            chic_taps (Data with digested by mnase enzyme and methylation converted by TAPS), chic_taps_transcriptome (Same as chic_taps, but includes annotations)
            chic_nla
            scartrace  (lineage tracing protocol)


        custom_flags(str): Arguments passed to the query name flagger, comma separated "MI,RX,bi,SM"

        ref(str) : Path to reference fasta file, autodected from bam header when not supplied

        umi_hamming_distance(int) : Max hamming distance on UMI's

        head (int) : Amount of molecules to process

        contig (str) : only process this contig

        region_start(int) : Zero based start coordinate of single region to process

        region_end(int) : Zero based end coordinate of single region to process, None: all contigs when contig is not set, complete contig when contig is set.

        alleles (str) : path to allele VCF

        allele_samples(str): Comma separated samples to extract from the VCF file. For example B6,SPRET

        unphased_alleles(str) : Path to VCF containing unphased germline SNPs

        mapfile (str) : 'Path to \*.safe.bgzf file, used to decide if molecules are uniquely mappable, generate one using createMapabilityIndex.py

        annotmethod (int) : Annotation resolving method. 0: molecule consensus aligned blocks. 1: per read per aligned base

        cluster (bool) : Run contigs in separate cluster jobs

        resolve_unproperly_paired_reads(bool) : When enabled bamtagmultiome will look through the complete bam file in a hunt for the mate, the two mates will always end up in 1 molecule if both present in the bam file. This also works when the is_proper_pair bit is not set. Use this option when you want to find the breakpoints of genomic re-arrangements.

        no_rejects(bool) : Do not write rejected reads

        mem (int) : Amount of gigabytes to request for cluster jobs

        time(int) : amount of wall clock hours to request for cluster jobs

        exons(str): Path to exon annotation GTF file

        introns(str): Path to intron annotation GTF file

        consensus(bool) : Calculate molecule consensus read, this feature is _VERY_ experimental

        consensus_model(str) : Path to consensus calling model, when none specified, this is learned based on the supplied bam file, ignoring sites supplied by -consensus_mask_variants

        consensus_mask_variants(str): Path VCF file masked for training on consensus caller

        consensus_n_train(int) : Amount of bases used for training the consensus model

        no_source_reads(bool) :  Do not write original reads, only consensus

        scartrace_r1_primers(str) : comma separated list of R1 primers used in scartrace protocol


    """

    MISC_ALT_CONTIGS_SCMO = 'MISC_ALT_CONTIGS_SCMO'
    every_fragment_as_molecule = args.every_fragment_as_molecule
    skip_contig = set(args.skip_contig.split(',')) if args.skip_contig is not None else set()


    if not args.o.endswith('.bam'):
        raise ValueError(
            "Supply an output which ends in .bam, for example -o output.bam")

    write_status(args.o,'unfinished')

    # verify wether the input is SAM, if so we need to convert:
    if args.bamin.endswith('.sam'):
        print('Input file is in SAM format. Performing conversion.')
        bam_path =  args.bamin.replace('.sam','.bam')
        args.bamin, _ = bam_path, sam_to_bam(args.bamin, bam_path)

    # Verify wether the input file is indexed and sorted...
    if not args.ignore_bam_issues:
        verify_and_fix_bam(args.bamin)

    for remove_existing_path in [args.o, f'{args.o}.bai']:
        if os.path.exists(remove_existing_path):
            print(f"Removing existing file {remove_existing_path}")
            os.remove(remove_existing_path)

    input_bam = pysam.AlignmentFile(args.bamin, "rb", ignore_truncation=args.ignore_bam_issues, threads=4)

    # autodetect reference:
    reference = None
    if args.ref is None:
        args.ref = get_reference_from_pysam_alignmentFile(input_bam)

    if args.ref is not None:
        try:
            reference = UnitialisedClass(CachedFastaNoHandle, args.ref)
            print(f'Loaded reference from {args.ref}')
        except Exception as e:
            print("Error when loading the reference file, continuing without a reference")
            reference = None

    ##### Set up consensus
    consensus_model_args = {'consensus_mode':None}
    if args.consensus:
        consensus_model_args = {'consensus_mode': 'majority'}
    if args.no_source_reads:
        consensus_model_args['no_source_reads'] = True

    ##### Define fragment and molecule class arguments and instances: ####

    query_name_flagger = None
    if args.qflagger is not None:
        if args.qflagger == 'custom_flags':
            query_name_flagger = CustomAssingmentQueryNameFlagger(
                args.custom_flags.split(','))
        elif args.qflagger == 'bulk':
            query_name_flagger = BulkFlagger()

        else:
            raise ValueError("Select from 'custom_flags, ..' ")


    molecule_class_args = {
        'umi_hamming_distance': args.umi_hamming_distance,
        'reference': reference
    }

    fragment_class_args = {
        'read_group_format' : args.read_group_format

    }
    yield_invalid = True  # if invalid reads should be written
    yield_overflow = True  # if overflow reads should be written

    if args.max_fragment_size is not None:
        fragment_class_args['max_fragment_size'] = args.max_fragment_size

    if args.no_rejects:
        yield_invalid = False

    if args.no_overflow:
        yield_overflow = False


    ignore_conversions = None
    if args.method == 'nla_taps' or args.method == 'chic_taps':
        ignore_conversions = set([('C', 'T'), ('G', 'A')])

    if args.alleles is not None and args.alleles!='none':
        molecule_class_args['allele_resolver'] = singlecellmultiomics.alleleTools.AlleleResolver(
            args.alleles,
            select_samples=args.allele_samples.split(',') if args.allele_samples is not None else None,
            lazyLoad=True,
            use_cache=args.use_allele_cache,
            verbose = args.set_allele_resolver_verbose,
            ignore_conversions=ignore_conversions)

    if args.mapfile is not None:
        assert args.mapfile.endswith('safe.bgzf'), 'the mapfile name should end with safe.bgzf '
        molecule_class_args['mapability_reader'] = MapabilityReader(args.mapfile)

    ### Transcriptome configuration ###
    if args.method in ('cs', 'vasa',  'chict') or args.method.endswith('_transcriptome'):
        print(
            colorama.Style.BRIGHT +
            'Running in transcriptome annotation mode' +
            colorama.Style.RESET_ALL)
        if args.exons is None :
            raise ValueError("Supply an exon GTF file")

        if args.introns is not None and args.exons is None:
            raise ValueError("Please supply both intron and exon GTF files")

        transcriptome_features = singlecellmultiomics.features.FeatureContainer(verbose=args.feature_container_verbose)
        print("Loading exons", end='\r')
        transcriptome_features.preload_GTF(
            path=args.exons,
            select_feature_type=['exon'],
            identifierFields=(
                'exon_id',
                'gene_id'),
            store_all=True,
            contig=args.contig,
            head=None)

        if args.introns is not None:
            print("Loading introns", end='\r')
            transcriptome_features.preload_GTF(
                path=args.introns,
                select_feature_type=['intron'],
                identifierFields=['transcript_id'],
                store_all=True,
                contig=args.contig,
                head=None)
        print("All features loaded")

        # Add more molecule class arguments
        transcriptome_feature_args = {
            'features': transcriptome_features,
            'auto_set_intron_exon_features': True,
        }


    bp_per_job = 10_000_000
    pooling_method=1
    bp_per_segment = 999_999_999 #@todo make this None or so
    fragment_size = 500
    one_contig_per_process=False
    max_time_per_segment = args.max_time_per_segment

    ### Method specific configuration ###
    if args.method == 'qflag':

        molecule_class = singlecellmultiomics.molecule.Molecule
        fragment_class = singlecellmultiomics.fragment.FragmentStartPosition

        # Write all reads
        yield_invalid = True

        bp_per_job = 3_000_000
        bp_per_segment = 3_000_000
        fragment_size = 0

    elif args.method == 'chic':
        molecule_class = singlecellmultiomics.molecule.CHICMolecule
        fragment_class = singlecellmultiomics.fragment.CHICFragment

        bp_per_job = 5_000_000
        bp_per_segment = 50_000
        fragment_size = 1000
        if max_time_per_segment is None:
            max_time_per_segment = 20*60 #



    elif args.method == 'nla' or args.method == 'nla_no_overhang':
        molecule_class = singlecellmultiomics.molecule.NlaIIIMolecule
        fragment_class = singlecellmultiomics.fragment.NlaIIIFragment

        if args.method == 'nla_no_overhang':
            assert reference is not None, 'Supply a reference fasta using -ref!'
            fragment_class_args.update({
                    'reference': reference,
                    'no_overhang': True
                })
        bp_per_job = 10_000_000
        bp_per_segment = 10_000_000
        fragment_size = 1000

    elif args.method == 'chic_nla':
        molecule_class=singlecellmultiomics.molecule.CHICNLAMolecule
        fragment_class=singlecellmultiomics.fragment.CHICFragment
        assert reference is not None, 'Supply a reference fasta using -ref!'
        molecule_class_args.update({
                'reference': reference,
        })

        bp_per_job = 5_000_000
        bp_per_segment = 50_000
        fragment_size = 1000

    elif args.method == 'cs_feature_counts' :
        molecule_class = singlecellmultiomics.molecule.Molecule
        fragment_class = singlecellmultiomics.fragment.FeatureCountsSingleEndFragment

    elif args.method == 'fl_feature_counts':

        molecule_class = singlecellmultiomics.molecule.Molecule
        fragment_class = singlecellmultiomics.fragment.FeatureCountsFullLengthFragment

    elif args.method == 'episeq':
        molecule_class = singlecellmultiomics.molecule.Molecule
        fragment_class = singlecellmultiomics.fragment.FeatureCountsSingleEndFragment

    elif args.method == 'nla_transcriptome':
        molecule_class = singlecellmultiomics.molecule.AnnotatedNLAIIIMolecule
        fragment_class = singlecellmultiomics.fragment.NlaIIIFragment

        molecule_class_args.update(transcriptome_feature_args)
        molecule_class_args.update({
            'stranded': None  # data is not stranded
        })
    elif args.method == 'chict' :

        molecule_class = singlecellmultiomics.molecule.AnnotatedCHICMolecule
        fragment_class = singlecellmultiomics.fragment.CHICFragment
        bp_per_job = 100_000_000 # One contig at a time to prevent reading the annotations too much
        bp_per_segment = 500_000
        fragment_size = 50_000

        molecule_class_args.update(transcriptome_feature_args)


    elif args.method == 'nla_taps':
        molecule_class = singlecellmultiomics.molecule.TAPSNlaIIIMolecule
        fragment_class = singlecellmultiomics.fragment.NlaIIIFragment

        molecule_class_args.update({
            'reference': reference,
            'taps': singlecellmultiomics.molecule.TAPS()
        })
        if args.consensus:
            bp_per_job = 2_000_000
        else:
            bp_per_job = 5_000_000
        bp_per_segment = 1_000_000
        fragment_size = 1000
    elif args.method == 'nla_taps_transcriptome': # Annotates reads in transcriptome
        molecule_class = singlecellmultiomics.molecule.AnnotatedTAPSNlaIIIMolecule
        fragment_class = singlecellmultiomics.fragment.NlaIIIFragment

        molecule_class_args.update(transcriptome_feature_args)

        molecule_class_args.update({
            'reference': reference,
            'taps': singlecellmultiomics.molecule.TAPS()
        })
        bp_per_job = 5_000_000
        bp_per_segment = 5_000_000
        fragment_size = 100_000

    elif args.method == 'nla_tapsp_transcriptome': # Annotates reads in transcriptome
        molecule_class = singlecellmultiomics.molecule.TAPSPTaggedMolecule
        fragment_class = singlecellmultiomics.fragment.NlaIIIFragment
        molecule_class_args.update(transcriptome_feature_args)
        molecule_class_args.update({
            'reference': reference,
            'taps': singlecellmultiomics.molecule.TAPS()
        })
        bp_per_job = 5_000_000
        bp_per_segment = 5_000_000
        fragment_size = 100_000

    elif args.method == 'chic_taps_transcriptome':

        bp_per_job = 5_000_000
        bp_per_segment = 5_000_000
        fragment_size = 100_000
        molecule_class_args.update({
            'reference': reference,
            'taps': singlecellmultiomics.molecule.TAPS(),
            'taps_strand':'R'
        })
        molecule_class_args.update(transcriptome_feature_args)
        molecule_class = singlecellmultiomics.molecule.AnnotatedTAPSCHICMolecule
        fragment_class = singlecellmultiomics.fragment.CHICFragment

    elif args.method == 'chic_taps':
        bp_per_job = 5_000_000
        bp_per_segment = 1_000_000
        fragment_size = 500
        molecule_class_args.update({
            'reference': reference,
            'taps': singlecellmultiomics.molecule.TAPS(),
            'taps_strand':'R'
        })
        molecule_class = singlecellmultiomics.molecule.TAPSCHICMolecule
        fragment_class = singlecellmultiomics.fragment.CHICFragment

    elif args.method == 'vasa' or args.method == 'cs':
        one_contig_per_process=True
        bp_per_job = 5_000_000
        bp_per_segment = 1_000_000
        fragment_size = 100_000

        molecule_class = singlecellmultiomics.molecule.TranscriptMolecule
        fragment_class = singlecellmultiomics.fragment.SingleEndTranscriptFragment
        #molecule_class_args.update(transcriptome_feature_args)
        fragment_class_args.update(transcriptome_feature_args)
        fragment_class_args.update({'stranded': True })
        molecule_class_args.update({'cache_size':1000})


    elif args.method == 'scartrace':

        molecule_class = singlecellmultiomics.molecule.ScarTraceMolecule
        fragment_class = singlecellmultiomics.fragment.ScarTraceFragment

        r1_primers = args.scartrace_r1_primers.split(',')
        fragment_class_args.update({
                'scartrace_r1_primers': r1_primers,
                #'reference': reference
            })


    else:
        raise ValueError("Supply a valid method")

    # Allow or disallow cycle shift:
    if args.allow_cycle_shift and fragment_class is singlecellmultiomics.fragment.NlaIIIFragment:
        fragment_class_args['allow_cycle_shift'] = True


    # This disables restriction motif checking
    if args.no_restriction_motif_check:
        fragment_class_args['check_motif'] = False

    # This disables umi_cigar_processing:
    if args.no_umi_cigar_processing:
        fragment_class_args['no_umi_cigar_processing'] = True

    if args.max_associated_fragments is not None:
        molecule_class_args['max_associated_fragments'] = args.max_associated_fragments


    if args.assignment_radius is not None:
        fragment_class_args['assignment_radius'] = args.assignment_radius
        if args.multiprocess:
            one_contig_per_process=True
            #raise NotImplementedError('-assignment_radius is currently incompatible with --multiprocess')

    if args.libname is not None:
        fragment_class_args['library_name'] = args.libname

    # This decides what molecules we will traverse
    if args.contig == MISC_ALT_CONTIGS_SCMO:
        contig = None
    else:
        contig = args.contig

    # This decides to only extract a single genomic region:
    if args.region_start is not None:
        if args.region_end is None:
            raise ValueError('When supplying -region_start then also supply -region_end')
        region_start = args.region_start
        region_end = args.region_end
    else:
        region_start = None
        region_end = None

    # Overwrite this flag when it is set:
    if args.one_contig_per_process:
        one_contig_per_process=True

    last_update = datetime.now()
    init_time = datetime.now()
    if args.molecule_iterator_verbosity_interval is not None and (args.molecule_iterator_verbose or (args.stats_file_path is not None )):

        stats_handle = None
        if args.stats_file_path is not None:
            stats_handle = open(args.stats_file_path,'w')

        def progress_callback_function( iteration, mol_iter, reads ):
            nonlocal last_update
            nonlocal init_time
            nonlocal stats_handle

            now = datetime.now()
            diff = (datetime.now()-last_update).total_seconds()
            if diff>args.molecule_iterator_verbosity_interval:

                diff_from_init = (datetime.now()-init_time).total_seconds()
                _contig, _pos = None, None
                for read in reads:
                    if read is not None:
                        _contig, _pos = read.reference_name, read.reference_start

                if args.molecule_iterator_verbose:
                    print( f'{mol_iter.yielded_fragments} fragments written, {mol_iter.deleted_fragments} fragments deleted ({(mol_iter.deleted_fragments/(mol_iter.deleted_fragments + mol_iter.yielded_fragments))*100:.2f} %), current pos: {_contig}, {_pos}, {mol_iter.waiting_fragments} fragments waiting             ' , end='\r')
                if stats_handle is not None:
                    stats_handle.write(f'{diff_from_init}\t{mol_iter.waiting_fragments}\t{mol_iter.yielded_fragments}\t{mol_iter.deleted_fragments}\t{_contig}\t{_pos}\n')
                    stats_handle.flush()
                last_update = now

    else:
        progress_callback_function = None


    molecule_iterator_args = {
        #'alignments': input_bam,
        'query_name_flagger': query_name_flagger,
        'molecule_class': molecule_class,
        'fragment_class': fragment_class,
        'molecule_class_args': molecule_class_args,
        'fragment_class_args': fragment_class_args,
        'yield_invalid': yield_invalid,
        'yield_overflow': yield_overflow,
        'start': region_start,
        'end': region_end,
        'contig': contig,
        'every_fragment_as_molecule': every_fragment_as_molecule,
        'skip_contigs':skip_contig,
        'progress_callback_function':progress_callback_function,
        'pooling_method' : pooling_method,
        'perform_allele_clustering': args.haplo_molecule_assignment and molecule_class_args.get('allele_resolver', None) is not None
    }



    if args.resolve_unproperly_paired_reads:
        molecule_iterator_args['iterator_class'] = MatePairIteratorIncludingNonProper

    if args.contig == MISC_ALT_CONTIGS_SCMO:
        # When MISC_ALT_CONTIGS_SCMO is set as argument, all molecules with reads
        # mapping to a contig returning True from the is_main_chromosome
        # function are used

        def Misc_contig_molecule_generator(input_bam, **molecule_iterator_args):
            for reference in input_bam.references:
                if not is_main_chromosome(reference):
                    molecule_iterator_args['contig'] = reference
                    yield from MoleculeIterator(input_bam, **molecule_iterator_args)

        molecule_iterator = Misc_contig_molecule_generator
    else:
        molecule_iterator = MoleculeIterator

    if args.method == 'qflag':
        molecule_iterator_args['every_fragment_as_molecule'] = True
        molecule_iterator_args['iterator_class'] = ReadIterator


    #####
    consensus_model_path = None


    # We needed to check if every argument is properly placed. If so; the jobs
    # can be sent to the cluster

    if args.cluster:
        if args.contig is None:
            write_status(args.o,'Submitting jobs. If this file remains, a job failed.')
            # Create jobs for all chromosomes:
            unique_id = str(uuid.uuid4())
            temp_prefix = os.path.abspath(os.path.dirname(
                args.o)) + '/SCMO_' + unique_id
            hold_merge = []

            ## Create folder to store cluster files:
            if args.clusterdir is None:
                cluster_file_folder = os.path.abspath(os.path.dirname(
                    args.o)) + '/cluster'
            else:
                cluster_file_folder = args.clusterdir
            print(f'Writing cluster scripts and standard out and error to {cluster_file_folder}')
            if not os.path.exists(cluster_file_folder):
                try:
                    os.makedirs(cluster_file_folder,exist_ok=True)
                except Exception as e:
                    print(e)
                    pass

            found_alts = 0
            files_to_merge = []
            for ci,chrom in enumerate([_chrom  for _chrom in
                        (list(input_bam.references) + [MISC_ALT_CONTIGS_SCMO])
                        if not _chrom in skip_contig]):

                if not is_main_chromosome(chrom):
                    found_alts += 1
                    continue
                if chrom == MISC_ALT_CONTIGS_SCMO and found_alts == 0:
                    continue

                temp_bam_path = f'{temp_prefix}_{chrom}.bam'

                if os.path.exists(temp_bam_path):
                    print(f"Removing existing temporary file {temp_bam_path}")
                    os.remove(temp_bam_path)

                arguments = " ".join(
                    [x for x in sys.argv if not x == args.o and x != '-o']) + f" -contig {chrom} -o {temp_bam_path}"
                files_to_merge.append(temp_bam_path)
                if consensus_model_path is not None:
                    arguments += f' -consensus_model {consensus_model_path}'
                job = f'SCMULTIOMICS_{ci}_{unique_id}'
                write_status(temp_bam_path,'SUBMITTED')
                job_id = submit_job(f'{arguments};', job_name=job, target_directory=cluster_file_folder,  working_directory=None,
                               threads_n=1, memory_gb=args.mem, time_h=args.time, scheduler=args.sched, copy_env=True,
                               email=None, mail_when_finished=False, hold=None,submit=True)


                print(f'Job for contig {chrom} submitted with job id: {job_id}')
                hold_merge.append(job_id)

            hold = hold_merge

            job = f'SCMULTIOMICS_MERGE_{unique_id}'

            if args.sched == 'local':
                hold = None

            final_status = args.o.replace('.bam','.status.txt')
            # Create list of output files
            command = f'samtools merge -@ 4 -c {args.o} {" ".join(files_to_merge)} && samtools index {args.o} && rm {temp_prefix}*.ba* && rm {temp_prefix}*.status.txt && echo "All done" > {final_status}'

            final_job_id = submit_job(f'{command};', job_name=job, target_directory=cluster_file_folder,  working_directory=None,
                           threads_n=4, memory_gb=10, time_h=args.time, scheduler=args.sched, copy_env=True,
                           email=None, mail_when_finished=False, hold=hold,submit=True)
            print(f'final job id is:{final_job_id}')
            exit()

    #####
    # Load unphased variants to memory
    """
    unphased_allele_resolver = None
    if args.unphased_alleles is not None:
        unphased_allele_resolver = singlecellmultiomics.alleleTools.AlleleResolver(
            use_cache=args.use_allele_cache,
            phased=False, ignore_conversions=ignore_conversions,verbose = args.set_allele_resolver_verbose)
        try:
            for i, variant in enumerate(
                pysam.VariantFile(
                    args.unphased_alleles).fetch(
                    args.contig)):
                if 'PASS' not in list(variant.filter):
                    continue
                if not all(
                        len(allele) == 1 for allele in variant.alleles) or len(
                        variant.alleles) != 2:
                    continue
                if sum([len(set(variant.samples[sample].alleles))
                        == 2 for sample in variant.samples]) < 2:
                    # Not heterozygous
                    continue

                unphased_allele_resolver.locationToAllele[variant.chrom][variant.pos - 1] = {
                    variant.alleles[0]: {'U'}, variant.alleles[1]: {'V'}}
        except Exception as e:  # todo catch this more nicely
            print(e)
    """



    if args.multiprocess:

        print("Tagging using multi-processing")
        tag_multiome_multi_processing(input_bam_path=args.bamin, out_bam_path=args.o, molecule_iterator=molecule_iterator,
                                      molecule_iterator_args=molecule_iterator_args,ignore_bam_issues=args.ignore_bam_issues,
                                      head=args.head, no_source_reads=args.no_source_reads,
                                      fragment_size=fragment_size, blacklist_path=args.blacklist,bp_per_job=bp_per_job,
                                      bp_per_segment=bp_per_segment, temp_folder_root=args.temp_folder, max_time_per_segment=max_time_per_segment,
                                      additional_args=consensus_model_args, n_threads=args.tagthreads, one_contig_per_process=one_contig_per_process,
                                      job_bed_file=args.jobbed
                                      )
    else:

        if consensus_model_args.get('consensus_mode') is not None:
            raise NotImplementedError('Please use --multiprocess')

        # Alignments are passed as pysam handle:
        if args.blacklist is not None:
            raise NotImplementedError("Blacklist can only be used with --multiprocess")
        tag_multiome_single_thread(
            args.bamin,
            args.o,
            molecule_iterator = molecule_iterator,
            molecule_iterator_args = molecule_iterator_args,
            consensus_model = None ,
            consensus_model_args={},
            ignore_bam_issues=False,
            head=args.head,
            no_source_reads=args.no_source_reads
            )


if __name__ == '__main__':
    args = argparser.parse_args()
    run_multiome_tagging(args)
