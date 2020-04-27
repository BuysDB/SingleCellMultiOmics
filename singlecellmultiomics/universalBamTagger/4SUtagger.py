#!/usr/bin/env python3
import pysam
from singlecellmultiomics.molecule import FourThiouridine, MoleculeIterator
from singlecellmultiomics.fragment import FeatureCountsSingleEndFragment

import singlecellmultiomics.pyutils
import pysamiterators
import collections
import glob
import pickle
import pandas as pd
import argparse

from singlecellmultiomics.bamProcessing.bamFunctions import sorted_bam_file, sort_and_index, get_reference_from_pysam_alignmentFile, add_readgroups_to_header, write_program_tag

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
    4FU detection tool
    """)
    argparser.add_argument('bamfile')
    argparser.add_argument('-head', type=int)
    argparser.add_argument(
        '-ref',
        type=str,
        help='reference (autodetected if not supplied)')
    args = argparser.parse_args()


    with pysam.AlignmentFile(args.bamfile, "rb") as alignments:

        reference = None
        if args.ref is None:
            args.ref = get_reference_from_pysam_alignmentFile(alignments)

        if args.ref is not None:
            reference = pysamiterators.iterators.CachedFasta(
                pysam.FastaFile(args.ref))

        fragment_class_args = {'umi_hamming_distance': 0}
        molecule_class_args = {'reference': reference}
        converted_base_frequencies = collections.defaultdict(
            lambda: collections.defaultdict(collections.Counter))

        with sorted_bam_file(f'{args.bamfile}.4SU_tagged.bam', origin_bam=alignments) as out:
            test_iter = MoleculeIterator(
                alignments,
                fragment_class_args=fragment_class_args,
                molecule_class_args={
                    'reference': reference},
                molecule_class=FourThiouridine,
                fragment_class=FeatureCountsSingleEndFragment)
            for i, molecule in enumerate(test_iter):
                molecule.write_tags()
                molecule.set_meta('mi', i)
                molecule.write_pysam(out)
                converted_base_frequencies[molecule.sample][molecule.gene][molecule.converted_bases] += 1
                if args.head is not None and i > (args.head - 1):
                    break

        for sample, data in converted_base_frequencies.items():
            pd.DataFrame(data).to_csv(f'{args.bamfile}.{sample}.counts.csv')
