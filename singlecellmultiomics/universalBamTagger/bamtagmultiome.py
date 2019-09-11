#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
from singlecellmultiomics.molecule import MoleculeIterator
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
from singlecellmultiomics.bamProcessing.bamFunctions import sort_and_index
from singlecellmultiomics.bamProcessing.bamFunctions import add_readgroups_to_header
import singlecellmultiomics.alleleTools
from singlecellmultiomics.universalBamTagger.customreads  import CustomAssingmentQueryNameFlagger

import argparse

argparser = argparse.ArgumentParser(
 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
 description='Assign molecules, set sample tags, set alleles')

argparser.add_argument('bamin',  type=str)
argparser.add_argument('-o',  type=str, help="output bam file", required=True)
argparser.add_argument('-qflagger',  type=str, default=None, help="Query flagging algorithm")
argparser.add_argument('-custom_flags',  type=str, default="MI,RX,BI,SM" )
argparser.add_argument('-method',  type=str, default=None, help="Protocol to tag, select from:nla,qflag,chic,nla_transcriptome,vasa,cs")
argparser.add_argument('-ref',  type=str, default=None, help="Path to reference fasta")
argparser.add_argument('-umi_hamming_distance',  type=int, default=1)
argparser.add_argument('-head',  type=int)
argparser.add_argument('-alleles',  type=str, help="Allele file (VCF)" )
argparser.add_argument('-allele_samples',  type=str, help="Comma separated samples to extract from the VCF file. For example B6,SPRET" )



argparser.add_argument('-annotmethod',  type=int, default=1, help="Annotation resolving method. 0: molecule consensus aligned blocks. 1: per read per aligned base" )
args = argparser.parse_args()
# autodetect reference:
reference = None
if args.ref is None:
    args.ref = get_reference_from_pysam_alignmentFile(input_bam)

if args.ref is not None:
    reference = pysamiterators.iterators.CachedFasta( pysam.FastaFile(args.ref) )


##### Define fragment and molecule class arguments and instances: ####

queryNameFlagger = None
if args.qflagger is not None:
    if args.qflagger == 'custom_flags':
        queryNameFlagger = CustomAssingmentQueryNameFlagger(args.custom_flags.split(','))
    else:
        raise ValueError("Select from 'custom_flags, ..' ")

molecule_class_args = {
    'umi_hamming_distance' : args.umi_hamming_distance,
    'reference' : args.ref
}
fragment_class_args = {}
yield_invalid= None # if invalid reads should be written

if args.alleles is not None:
    molecule_class_args['allele_resolver'] = \
        singlecellmultiomics.alleleTools.AlleleResolver(args.alleles,
                                                select_samples=args.allele_samples.split(',') if args.allele_samples is not None else None,
                                                lazyLoad=True )

if args.method=='qflag':
    moleculeClass = singlecellmultiomics.molecule.Molecule
    fragmentClass = singlecellmultiomics.fragment.Fragment
    # Write all reads
    yield_invalid = True

elif args.method=='chic':
    moleculeClass = singlecellmultiomics.molecule.CHICMolecule
    fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment

elif args.method=='nla':
    moleculeClass = singlecellmultiomics.molecule.NlaIIIMolecule
    fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment

elif args.method=='nla_transcriptome':
    moleculeClass = singlecellmultiomics.molecule.AnnotatedNLAIIIMolecule
    fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment

    molecule_class_args.update({
        'pooling_method' : 1, # all data from the same cell can be dealt with separately
        'stranded': None # data is not stranded
    })

elif args.method=='vasa' or args.method=='cs':
    moleculeClass = singlecellmultiomics.molecule.VASA
    fragmentClass = singlecellmultiomics.fragment.SingleEndTranscript

    molecule_class_args.update({
        'pooling_method' : 0, # all data from the same cell can be dealt with separately
        'stranded': 1 # data is not stranded
    })

else:
    raise ValueError("Supply a valid method")

# We needed to check if every argument is properly placed. If so; the jobs can be sent to the cluster
if args.cluster:
    if args.contig is None:
        # Create jobs for all chromosomes:
        temp_prefix = os.path.abspath( os.path.dirname(args.o) )+ '/' + str(uuid.uuid4())
        hold_merge=[]
        for chrom in input_bam.references:
            if chrom.startswith('KN') or chrom.startswith('KZ') or chrom.startswith('chrUn') or chrom.endswith('_random') or 'ERCC' in chrom  or chrom.endswith('_alt') or "HLA-" in chrom:
                continue
            temp_bam_path = f'{temp_prefix}_{chrom}.bam'
            arguments = " ".join([x for x in sys.argv if not x==args.o in x and x!='-o'])  + f" -contig {chrom} -o {temp_bam_path}"
            job = f'TAPS_{str(uuid.uuid4())}'
            os.system( f'submission.py --silent' + f' -y --py36 -time {args.time} -t 1 -m {args.mem} -N {job} " {arguments};"' )
            hold_merge.append(job)

        hold =  ','.join(hold_merge)
        os.system( f'submission.py --silent' + f' -y --py36 -time {args.time} -t 1 -m 10 -N {job} -hold {hold} " samtools merge {args.o} {temp_prefix}*.bam; samtools index {args.o}; rm {temp_prefix}*.ba*"' )
        exit()

#####

out_bam_path = args.o

# Temp bam file to write tagged records to. This file does not have read groups yet,
# is unsorted and has the same header as the input bam file
out_bam_temp_path = f'{out_bam_path}.unsorted'

# This is the path with read groups added:
out_bam_temp_path_rg = f'{out_bam_path}.unsorted.rg'

# Open the input bam file
with pysam.AlignmentFile(input_bam_path, "rb") as input_bam:
    #Copy the header
    input_header = input_bam.header.copy()

    # No "with" statement , because the nesting is _really_ ugly.
    # Sorry. Dont forget to close this handle. See: @close
    out_bam_temp = pysam.AlignmentFile(out_bam_temp_path, "wb", header = input_header)

    read_groups = set() # Store unique read groups in this set
    for i,molecule in enumerate(
            MoleculeIterator(
                alignments=input_bam,
                queryNameFlagger=queryNameFlagger,
                moleculeClass=moleculeClass,
                fragmentClass=fragmentClass,
                molecule_class_args=molecule_class_args,
                fragment_class_args=fragment_class_args,
                yield_invalid=yield_invalid
            )
        ):

        # Stop when enough molecules are processed
        if args.head is not None and (i-1)>=args.head:
            break

        # set unique molecule identifier
        molecule.set_meta('mi', i)
        # Write tag values
        molecule.write_tags()

        # Update read groups
        for fragment in molecule:
            read_groups.add(fragment.get_read_group())

        # Write the reads to the output file
        molecule.write_pysam( out_bam_temp )

# Add readgroups to the bam file
add_readgroups_to_header(
    out_bam_temp_path,
    read_groups,
    target_bam_path=out_bam_temp_path )

# Sort and index
sort_and_index( out_bam_temp_path,  out_bam_path, remove_unsorted=True)
