#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pysam
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import singlecellmultiomics.molecule
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
            try:
                call_dict = taps.molecule_to_context_call_dict(molecule)
            except ValueError as e:
                # Safe calls cannot be made! (R2 didn't map properly for any of the fragments)
                continue
            # Write bismark tags:
            molecule.set_methylation_call_tags(call_dict)

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
