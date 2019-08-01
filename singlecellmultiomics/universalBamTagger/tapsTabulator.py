#!/usr/bin/env python
# -*- coding: utf-8 -*-
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysamiterators
import pysam
import argparse
import singlecellmultiomics.bamProcessing.bamFunctions as bf

if __name__=='__main__':
    argparser = argparse.ArgumentParser(
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     description='Extract TAPS methylation calls from BAM file')

    argparser.add_argument('alignmentfile',  type=str)
    argparser.add_argument('-ref',  type=str, help='path to reference fasta file, auto detected from bamfile')
    argparser.add_argument('-head',  type=int)
    argparser.add_argument('-contig',  type=str,help='contig to run on, all when not specified')
    argparser.add_argument('-moleculeNameSep',  type=str,help='Separator to use in molecule name', default=':')
    args = argparser.parse_args()
    alignments = pysam.AlignmentFile(args.alignmentfile)

    if args.ref is None:
        args.ref = bf.get_reference_from_pysam_alignmentFile(alignments)
        if args.ref is None:
            raise ValueError("Supply reference, -ref")

    reference = pysamiterators.iterators.CachedFasta( pysam.FastaFile(args.ref) )
    taps = singlecellmultiomics.molecule.TAPS(reference=reference)

    for i,molecule in  enumerate( singlecellmultiomics.molecule.MoleculeIterator(
        alignments=alignments,
        moleculeClass=singlecellmultiomics.molecule.TAPSMolecule,
        fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
        fragment_class_args={'umi_hamming_distance':1},
        molecule_class_args={'reference':reference,'taps':taps},
        contig=args.contig)):
        if args.head and i>=args.head:
            break
        # If calls cannot be obtained skip the molecule
        if molecule.methylation_call_dict is None:
            continue

        for (chromosome, location),call in molecule.methylation_call_dict.items():
            if call=='.': # Only print calls concerning C's
                continue
            print(f'{molecule.sample}{args.moleculeNameSep}{molecule.umi}{args.moleculeNameSep}{molecule.get_strand_repr()}\t{chromosome}\t{location}\t{call}')
