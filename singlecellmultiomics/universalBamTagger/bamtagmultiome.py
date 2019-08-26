#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam

from singlecellmultiomics.molecule import MoleculeIterator
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment

from singlecellmultiomics.bamProcessing.bamFunctions import sort_and_index
from singlecellmultiomics.bamProcessing.bamFunctions import add_readgroups_to_header

from types import SimpleNamespace
args = SimpleNamespace(
            o='./tagging_test.bam',
            bamin = './tagged_transcriptome.bam',
            head=100
)


input_bam_path = args.bamin
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
                alignments=input_bam
            )
        ):

        # Stop when enough molecules are processed
        if args.head is not None and (i-1)>=args.head:
            break

        # Write tag values
        molecule.write_tags()

        # Update read groups
        for fragment in molecule:
            read_groups.add(fragment.get_read_group())

        # Write the reads to the output file
        molecule.write_pysam( out_bam_temp )


    # @close
    out_bam_temp.close()

    # Add readgroups to the bam file
    add_readgroups_to_header(
        out_bam_temp_path,
        read_groups,
        target_bam_path=out_bam_temp_path )

    # Sort and index
    sort_and_index( out_bam_temp_path,  out_bam_path, remove_unsorted=True)
