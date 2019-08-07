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
    argparser.add_argument('-minmq',  type=int, default=50)
    argparser.add_argument('-contig',  type=str,help='contig to run on, all when not specified')
    argparser.add_argument('-moleculeNameSep',  type=str,help='Separator to use in molecule name', default=':')
    argparser.add_argument('-samples',  type=str,help='Samples to select, separate with comma. For example CellA,CellC,CellZ', default=None)
    argparser.add_argument('-context',  type=str,help='Contexts to select, separate with comma. For example Z,H,X', default=None)
    args = argparser.parse_args()
    alignments = pysam.AlignmentFile(args.alignmentfile)

    samples = None if args.samples is None else set(args.samples.split(','))
    contexts = None if args.context is None else set(
        [x.upper() for x in  args.context.split(',')] +
        [x.lower() for x in  args.context.split(',')])

    if args.ref is None:
        args.ref = bf.get_reference_from_pysam_alignmentFile(alignments)
        if args.ref is None:
            raise ValueError("Supply reference, -ref")

    reference = pysamiterators.iterators.CachedFasta( pysam.FastaFile(args.ref) )
    taps = singlecellmultiomics.molecule.TAPS(reference=reference)

    try:
        for i,molecule in  enumerate( singlecellmultiomics.molecule.MoleculeIterator(
            alignments=alignments,
            moleculeClass=singlecellmultiomics.molecule.TAPSNlaIIIMolecule,
            fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
            fragment_class_args={'umi_hamming_distance':1},
            molecule_class_args={
                'reference':reference,
                'site_has_to_be_mapped':True,
                'taps':taps
            },
            contig=args.contig)):

            if args.head and i>=args.head:
                break
                
            if not molecule.is_valid() or molecule.is_multimapped() or molecule.get_mean_mapping_qual()<args.minmq:
                continue

            # If calls cannot be obtained skip the molecule
            if molecule.methylation_call_dict is None:
                continue

            # Skip sample if not selected
            if samples is not None and molecule.sample not in samples:
                continue

            for (chromosome, location),call in molecule.methylation_call_dict.items():
                if call=='.': # Only print calls concerning C's
                    continue

                # Skip non-selected contexts
                if contexts is not None and call not in contexts:
                    continue

                print(f'{molecule.sample}{args.moleculeNameSep}{i}{args.moleculeNameSep}{molecule.umi}{args.moleculeNameSep}{molecule.get_strand_repr()}\t{chromosome}\t{location}\t{call}')
    except (KeyboardInterrupt,BrokenPipeError) as e:
        pass
