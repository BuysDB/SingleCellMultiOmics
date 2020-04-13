#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
from singlecellmultiomics.molecule import MoleculeIterator
import singlecellmultiomics
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
from singlecellmultiomics.bamProcessing.bamFunctions import sorted_bam_file, get_reference_from_pysam_alignmentFile, write_program_tag, MapabilityReader, verify_and_fix_bam

from singlecellmultiomics.utils import is_main_chromosome
from singlecellmultiomics.utils.submission import submit_job
import singlecellmultiomics.alleleTools
from singlecellmultiomics.universalBamTagger.customreads import CustomAssingmentQueryNameFlagger
import singlecellmultiomics.features
from pysamiterators import CachedFasta,MatePairIteratorIncludingNonProper,MatePairIterator

import argparse
import uuid
import os
import sys
import colorama
import sklearn
import pkg_resources
import pickle
from datetime import datetime
import traceback

available_consensus_models = pkg_resources.resource_listdir('singlecellmultiomics','molecule/consensus_model')

argparser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Assign molecules, set sample tags, set alleles')
argparser.add_argument('bamin', type=str)
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
    nla_no_overhang (Data with digested by Nla III enzyme, without the CATG present in the reads)
    scartrace (Lineage tracing )
    """)
argparser.add_argument(
    '-qflagger',
    type=str,
    default=None,
    help="Query flagging algorithm")
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
    '--use_allele_cache',
    action='store_true',
    help='Write and use a cache file for the allele information. NOTE: THIS IS NOT THREAD SAFE! Meaning you should not use this function on multiple libraries at the same time when the cache files are not available. Once they are available there is not thread safety issue anymore')

argparser.add_argument('-molecule_iterator_verbosity_interval',type=int,default=None,help='Molecule iterator information interval in seconds')
argparser.add_argument('--molecule_iterator_verbose', action='store_true', help='Show progress indication on command line')
argparser.add_argument('-stats_file_path',type=str,default=None,help='Path to logging file, ends with ".tsv"')

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
    '--no_rejects',
    action='store_true',
    help='Do not write rejected reads to output file')
fragment_settings.add_argument(
    '--no_overflow',
    action='store_true',
    help='Do not write overflow reads to output file. Overflow reads are reads which are discarded because the molecule reached the maximum capacity of associated fragments')


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
cg.add_argument('-consensus_mask_variants', type=str,
                help='variants to mask during training (VCF file)')
cg.add_argument(
    '-consensus_model',
    type=str,
    help=f'Name of or path to consensus classifier, built-in available models are:{", ".join(available_consensus_models)}',
    default=None)
cg.add_argument(
    '-consensus_n_train',
    type=int,
    help='Amount of bases used for training',
    default=500_000)

cg.add_argument(
    '-consensus_k_rad',
    type=int,
    help='consensus model k radius',
    default=3)


cg.add_argument('--no_source_reads', action='store_true',
                help='Do not write original reads, only consensus ')

cg.add_argument('--consensus_allow_train_location_oversampling', action='store_true',
                help='Allow to train the consensus model multiple times for a single genomic location')


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

        bamin (str) : bam file to process

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
            chic_taps (Data with digested by mnase enzyme and methylation converted by TAPS)
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
            reference = CachedFasta(
                pysam.FastaFile(args.ref))
            print(f'Loaded reference from {args.ref}')
        except Exception as e:
            print("Error when loading the reference file, continuing without a reference")
            reference = None

    ##### Define fragment and molecule class arguments and instances: ####

    queryNameFlagger = None
    if args.qflagger is not None:
        if args.qflagger == 'custom_flags':
            queryNameFlagger = CustomAssingmentQueryNameFlagger(
                args.custom_flags.split(','))
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

    if args.alleles is not None:
        molecule_class_args['allele_resolver'] = singlecellmultiomics.alleleTools.AlleleResolver(
            args.alleles,
            select_samples=args.allele_samples.split(',') if args.allele_samples is not None else None,
            lazyLoad=True,
            use_cache=args.use_allele_cache,
            verbose = args.set_allele_resolver_verbose,
            ignore_conversions=ignore_conversions)

    if args.mapfile is not None:
        molecule_class_args['mapability_reader'] = MapabilityReader(
            args.mapfile)

    ### Transcriptome configuration ###
    if args.method in ('nla_transcriptome', 'cs', 'vasa'):
        print(
            colorama.Style.BRIGHT +
            'Running in transcriptome annotation mode' +
            colorama.Style.RESET_ALL)
        if args.exons is None :
            raise ValueError("Supply an exon GTF file")

        if args.introns is not None and args.exons is None:
            raise ValueError("Please supply both intron and exon GTF files")

        transcriptome_features = singlecellmultiomics.features.FeatureContainer()
        print("Loading exons", end='\r')
        transcriptome_features.loadGTF(
            args.exons,
            select_feature_type=['exon'],
            identifierFields=(
                'exon_id',
                'gene_id'),
            store_all=True,
            contig=args.contig,
            head=None)

        if args.introns is not None:
            print("Loading introns", end='\r')
            transcriptome_features.loadGTF(
                args.introns,
                select_feature_type=['intron'],
                identifierFields=['transcript_id'],
                store_all=True,
                contig=args.contig,
                head=None)
        print("All features loaded")

        # Add more molecule class arguments
        molecule_class_args.update({
            'features': transcriptome_features,
            'auto_set_intron_exon_features': True
        })

    ### Method specific configuration ###
    if args.method == 'qflag':
        moleculeClass = singlecellmultiomics.molecule.Molecule
        fragmentClass = singlecellmultiomics.fragment.Fragment
        # Write all reads
        yield_invalid = True

    elif args.method == 'chic':
        moleculeClass = singlecellmultiomics.molecule.CHICMolecule
        fragmentClass = singlecellmultiomics.fragment.CHICFragment

    elif args.method == 'nla' or args.method == 'nla_no_overhang':
        moleculeClass = singlecellmultiomics.molecule.NlaIIIMolecule
        fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment

        if args.method == 'nla_no_overhang':
            assert reference is not None, 'Supply a reference fasta using -ref!'
            fragment_class_args.update({
                    'reference': reference,
                    'no_overhang': True
                })

    elif args.method == 'chic_nla':
        moleculeClass=singlecellmultiomics.molecule.CHICNLAMolecule
        fragmentClass=singlecellmultiomics.fragment.CHICFragment
        assert reference is not None, 'Supply a reference fasta using -ref!'
        molecule_class_args.update({
                'reference': reference,
        })

    elif args.method == 'cs_feature_counts' :
        moleculeClass = singlecellmultiomics.molecule.Molecule
        fragmentClass = singlecellmultiomics.fragment.FeatureCountsSingleEndFragment

    elif args.method == 'fl_feature_counts':

        moleculeClass = singlecellmultiomics.molecule.Molecule
        fragmentClass = singlecellmultiomics.fragment.FeatureCountsFullLengthFragment

    elif args.method == 'episeq' :
        moleculeClass = singlecellmultiomics.molecule.Molecule
        fragmentClass = singlecellmultiomics.fragment.FeatureCountsSingleEndFragment

    elif args.method == 'nla_transcriptome':
        moleculeClass = singlecellmultiomics.molecule.AnnotatedNLAIIIMolecule
        fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment

        molecule_class_args.update({
            'pooling_method': 1,  # all data from the same cell can be dealt with separately
            'stranded': None  # data is not stranded
        })

    elif args.method == 'nla_taps':
        moleculeClass = singlecellmultiomics.molecule.TAPSNlaIIIMolecule
        fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment

        molecule_class_args.update({
            'reference': reference,
            'taps': singlecellmultiomics.molecule.TAPS(reference=reference)
        })

    elif args.method == 'chic_taps':

        molecule_class_args.update({
            'reference': reference,
            'taps': singlecellmultiomics.molecule.TAPS(reference=reference)
        })
        moleculeClass = singlecellmultiomics.molecule.TAPSCHICMolecule
        fragmentClass = singlecellmultiomics.fragment.CHICFragment

    elif args.method == 'vasa' or args.method == 'cs':
        moleculeClass = singlecellmultiomics.molecule.VASA
        fragmentClass = singlecellmultiomics.fragment.SingleEndTranscript

        molecule_class_args.update({
            'pooling_method': 1,  # all data from the same cell can be dealt with separately
            'stranded': 1  # data is stranded
        })

    elif args.method == 'scartrace':

        moleculeClass = singlecellmultiomics.molecule.ScarTraceMolecule
        fragmentClass = singlecellmultiomics.fragment.ScarTraceFragment

        r1_primers = args.scartrace_r1_primers.split(',')
        fragment_class_args.update({
                'scartrace_r1_primers': r1_primers,
                #'reference': reference
            })


    else:
        raise ValueError("Supply a valid method")

    # Allow or disallow cycle shift:
    if args.allow_cycle_shift and fragmentClass is singlecellmultiomics.fragment.NLAIIIFragment:
        fragment_class_args['allow_cycle_shift'] = True

    # This disables umi_cigar_processing:
    if args.no_umi_cigar_processing:
        fragment_class_args['no_umi_cigar_processing'] = True

    if args.max_associated_fragments is not None:
        molecule_class_args['max_associated_fragments'] = args.max_associated_fragments

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
        'alignments': input_bam,
        'queryNameFlagger': queryNameFlagger,
        'moleculeClass': moleculeClass,
        'fragmentClass': fragmentClass,
        'molecule_class_args': molecule_class_args,
        'fragment_class_args': fragment_class_args,
        'yield_invalid': yield_invalid,
        'yield_overflow': yield_overflow,
        'start':region_start,
        'end':region_end,
        'contig': contig,
        'every_fragment_as_molecule': every_fragment_as_molecule,
        'skip_contigs':skip_contig,
        'progress_callback_function':progress_callback_function
    }

    if args.resolve_unproperly_paired_reads:
        molecule_iterator_args['iterator_class'] = MatePairIteratorIncludingNonProper

    if args.contig == MISC_ALT_CONTIGS_SCMO:
        # When MISC_ALT_CONTIGS_SCMO is set as argument, all molecules with reads
        # mapping to a contig returning True from the is_main_chromosome
        # function are used

        def Misc_contig_molecule_generator(molecule_iterator_args):
            for reference in input_bam.references:
                if not is_main_chromosome(reference):
                    molecule_iterator_args['contig'] = reference
                    yield from MoleculeIterator(**molecule_iterator_args)

        molecule_iterator = Misc_contig_molecule_generator(
            molecule_iterator_args)
    else:
        molecule_iterator = MoleculeIterator(**molecule_iterator_args)

    #####
    consensus_model_path = None

    if args.consensus:
        # Load from path if available:

        if args.consensus_model is not None:
            if os.path.exists(args.consensus_model):
                model_path = args.consensus_model
            else:
                model_path = pkg_resources.resource_filename(
                    'singlecellmultiomics', f'molecule/consensus_model/{args.consensus_model}')

            if model_path.endswith('.h5'):
                try:
                    from tensorflow.keras.models import load_model
                except ImportError:
                    print("Please install tensorflow")
                    raise
                consensus_model = load_model(model_path)

            else:
                with open(model_path, 'rb') as f:
                    consensus_model = pickle.load(f)
        else:
            skip_already_covered_bases = not args.consensus_allow_train_location_oversampling
            if args.consensus_mask_variants is None:
                mask_variants = None
            else:
                mask_variants = pysam.VariantFile(args.consensus_mask_variants)
            print("Fitting consensus model, this may take a long time")
            consensus_model = singlecellmultiomics.molecule.train_consensus_model(
                molecule_iterator,
                mask_variants=mask_variants,
                n_train=args.consensus_n_train,
                skip_already_covered_bases=skip_already_covered_bases
                )
            # Write the consensus model to disk
            consensus_model_path = os.path.abspath(
                os.path.dirname(args.o)) + '/consensus_model.pickle.gz'
            print(f'Writing consensus model to {consensus_model_path}')
            with open(consensus_model_path, 'wb') as f:
                pickle.dump(consensus_model, f)

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
    out_bam_path = args.o

    # Copy the header
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

    read_groups = dict()  # Store unique read groups in this dict
    with sorted_bam_file(out_bam_path, header=input_header, read_groups=read_groups) as out:

        try:
            for i, molecule in enumerate(molecule_iterator):

                # Stop when enough molecules are processed
                if args.head is not None and (i - 1) >= args.head:
                    break

                # set unique molecule identifier
                molecule.set_meta('mi', f'{molecule.get_a_reference_id()}_{i}')

                # Write tag values
                molecule.write_tags()

                if unphased_allele_resolver is not None:  # write unphased allele tag:
                    molecule.write_allele_phasing_information_tag(
                        unphased_allele_resolver, 'ua')

                # Update read groups
                for fragment in molecule:
                    rgid = fragment.get_read_group()
                    if not rgid in read_groups:
                        read_groups[rgid] = fragment.get_read_group(True)[1]

                # Calculate molecule consensus
                if args.consensus:
                    try:
                        consensus_reads = molecule.deduplicate_to_single_CIGAR_spaced(
                            out,
                            f'consensus_{molecule.get_a_reference_id()}_{i}',
                            consensus_model,
                            NUC_RADIUS=args.consensus_k_rad
                            )
                        for consensus_read in consensus_reads:
                            consensus_read.set_tag('RG', molecule[0].get_read_group())
                            consensus_read.set_tag('mi', i)
                            out.write(consensus_read)
                    except Exception as e:

                        #traceback.print_exc()
                        #print(e)
                        molecule.set_rejection_reason('CONSENSUS_FAILED',set_qcfail=True)
                        molecule.write_pysam(out)


                # Write the reads to the output file
                if not args.no_source_reads:
                    molecule.write_pysam(out)
        except Exception as e:
            write_status(args.o,'FAIL, The file is not complete')
            raise e

        # Reached the end of the generator
        write_status(args.o,'Reached end. All ok!')


if __name__ == '__main__':
    args = argparser.parse_args()
    run_multiome_tagging(args)
