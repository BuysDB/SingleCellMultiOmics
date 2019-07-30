#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pysam
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import gzip
import collections
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import pysamiterators
import sys
import os
import uuid
if __name__=='__main__':
    argparser = argparse.ArgumentParser(
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     description='Add methylation information to BAM file')

    argparser.add_argument('alignmentfile',  type=str)
    argparser.add_argument('-o',  type=str, help="output BAM path", required=True)

    argparser.add_argument('-ref',  type=str, required=True, help='path to reference fasta file ')

    argparser.add_argument('-head',  type=int)
    argparser.add_argument('-contig',  type=str,help='contig to run on')
    args = argparser.parse_args()


    g= pysam.AlignmentFile(args.alignmentfile)
    if args.contig is None:
        # Create jobs for all chromosomes:
        temp_prefix = os.path.abspath( os.path.dirname(args.o) )+ '/' + str(uuid.uuid4())
        hold_merge=[]
        for chrom in g.references:
            if chrom.startswith('chrUn') or chrom.endswith('_random'):
                continue
            temp_bam_path = f'{temp_prefix}_{chrom}.bam'
            arguments = " ".join([x for x in sys.argv if not x==args.o in x and x!='-o'])  + f" -contig {chrom} -o {temp_bam_path}"
            job = f'TAPS_{str(uuid.uuid4())}'
            print( f'submission.py' + f' -y --py36 -time 50 -t 1 -m 8 -N {job} " {arguments};"' )
            hold_merge.append(job)

        hold =  ','.join(hold_merge)
        print( f'submission.py' + f' -y --py36 -time 50 -t 1 -m 8 -N {job} -hold {hold} " samtools merge {args.o} {temp_prefix}*.bam; samtools index {args.o}; rm {temp_prefix}*.bam"' )
        exit()

    reference = pysamiterators.iterators.CachedFasta( pysam.FastaFile(args.ref) )
    taps = singlecellmultiomics.molecule.TAPS(reference=reference)
    temp_out = f'{args.o}.temp.out.bam'

    with pysam.AlignmentFile(temp_out , "wb",header=g.header) as output:
        for i,molecule in  enumerate( singlecellmultiomics.molecule.MoleculeIterator(
            alignments=g,
            moleculeClass=singlecellmultiomics.molecule.Molecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
            fragment_class_args={'umi_hamming_distance':1},
            contig=args.contig)):
            if args.head and i>args.head:
                break
            # Obtain taps methylation calls:

            if molecule.is_multimapped():
                continue

            # Find all aligned positions and corresponding reference bases:
            aligned_reference_positions = {} #(chrom,pos)->base
            for read in molecule.iter_reads():
                for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True, matches_only=True):
                    aligned_reference_positions[(read.reference_name,ref_pos)] = ref_base.upper()

            # Obtain consensus:
            try:
                consensus = molecule.get_consensus()
            except ValueError:
                continue # we cannot obtain a consensus


            # if insecure about the consensus emit the reference, if no CG AT conversion skip
            for location, reference_base in aligned_reference_positions.items():
                if not location in consensus or reference_base not in 'CG' or consensus[location] not in 'AT':
                    consensus[location] = reference_base.upper()

            # find all locations where a C/G was converted into A/T, now strand specific
            converted_bases = 0
            conversions = {}
            for location, reference_base in aligned_reference_positions.items():
                if (not molecule.strand and reference_base=='C' and consensus[location] in 'CT') or \
                    molecule.strand and reference_base=='G' and consensus[location] in 'AG':
                    conversions[location] = {'ref':reference_base, 'obs':consensus[location]}
                    if consensus[location] in 'TA':
                        converted_bases+=1

            # obtain the context of the conversions:
            conversion_contexts  ={
                location : taps.position_to_context(
                    *location,
                    observed_base = observations['obs'],
                    strand = molecule.strand)[1]
                for location, observations in conversions.items()}

            # Write bismark tags:
            molecule.set_methylation_call_tags(conversion_contexts)

            for read in molecule.iter_reads():
                if molecule.strand == 0: # forward
                    read.set_tag('XR','GA')
                    read.set_tag('XG','GA')
                else:
                    read.set_tag('XR','GA')
                    read.set_tag('XG','CT')
            
            # Write all reads to the new bam file:
            for read in molecule.iter_reads():
                if read is not None:
                    output.write(read)

    # Sort and index
    # Perform a reheading, sort and index
    cmd = f"""samtools sort {temp_out} > {args.o}; samtools index {args.o};
    rm {temp_out};
    """
    os.system(cmd)
