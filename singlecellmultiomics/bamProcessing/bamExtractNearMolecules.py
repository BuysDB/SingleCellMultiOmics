#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pysam
from array import array
from singlecellmultiomics.bamProcessing.bamAnalyzeCutDistances import get_stranded_pairwise_counts, get_sc_cut_dictionary, strict_read_counts_function
from singlecellmultiomics.bamProcessing.bamFunctions import get_contigs_with_reads, sorted_bam_file
from multiprocessing import Pool
from singlecellmultiomics.utils import is_main_chromosome
from collections import defaultdict
from pysamiterators import MatePairIterator
from singlecellmultiomics.molecule import MoleculeIterator, CHICMolecule
from singlecellmultiomics.fragment import CHICFragment
from singlecellmultiomics.bamProcessing.bamFunctions import get_reference_path_from_bam

def get_nearby_reads():

    pairing_indices = {}
    current_pair_id = 0

    with pysam.AlignmentFile(bam_path,threads=8) as alignments:
        for contig in contigs:

            for R1, R2 in MatePairIterator(alignments, contig=contig):
                if R1 is None or R1.is_qcfail or R1.mapping_quality<50:
                    continue

                # Obtain the cut site and strand
                cut_location = R1.get_tag('DS')
                cell = R1.get_tag( 'SM')
                strand = R1.is_reverse

                if not cell in selection:
                    continue

                if not (R1.reference_name, cut_location) in selection[cell] :
                    continue

                # Obtain paired cut site
                paired_cut = pairings[(cell,contig,cut_location)]
                R1.set_tag('PC',int(paired_cut))
                if R2 is not None:
                    R2.set_tag('PC',int(paired_cut))

                # Set pair ID:
                k = (cell, contig, min(paired_cut, cut_location))
                if not k in pairing_indices:
                    pairing_indices[k] = current_pair_id
                    pid = current_pair_id
                    current_pair_id+=1
                else:
                    pid = pairing_indices.get(k,-1)

                R1.set_tag('BX',pid)
                if R2 is not None:
                    R2.set_tag('BX',pid)

                yield R1, R2


def crc(reads, molecule):
    bx = molecule[0][0].get_tag('BX')
    for r in reads:
        r.set_tag('BX',bx)


def extract_near_cuts(bam_path, output_path):

    cuts_stranded = get_sc_cut_dictionary(bam_path,strand_specific=True, filter_function=strict_read_counts_function)
    contigs = [contig for contig in get_contigs_with_reads(bam_path) if is_main_chromosome(contig) and contig[0]!='J' and contig[0]!='M']



if __name__ == '__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract nearby molecules. Sets BX tag for nearby molecules and PC (paired cut) location')

    argparser.add_argument('alignmentfile', type=str)
    argparser.add_argument('outputfile', type=str, help='output bam')

    args = argparser.parse_args()

    #extract_near_cuts(args.alignmentfile, args.outputfile)
    bam_path = args.alignmentfile
    output_path = args.outputfile
    cuts_stranded = get_sc_cut_dictionary(bam_path,strand_specific=True, filter_function=strict_read_counts_function)
    contigs = [contig for contig in get_contigs_with_reads(bam_path) if is_main_chromosome(contig) and contig[0]!='J' and contig[0]!='M']

    def create_molecule_selection(contig) :

        selection = defaultdict( set )
        pairings =  dict()


        for cell in cuts_stranded[contig]:
            fwd_cuts = np.array( sorted([pos for strand, pos in cuts_stranded[contig][cell] if not strand]) )
            rev_cuts = np.array( sorted([pos for strand, pos in cuts_stranded[contig][cell] if strand]) )


            for cut in fwd_cuts:
                distances = rev_cuts - cut

                near_idcs = (distances>0) & (distances<1500)

                if sum(near_idcs)>0:
                    for c in rev_cuts[near_idcs]:
                        selection[cell].add( (contig, c ))
                        pairings[(cell,contig,c)] = cut
                    selection[cell].add( (contig, cut ))
                    pairings[(cell,contig,cut)] = c


                continue
        return selection, pairings

    selection = defaultdict( set )
    pairings =  dict()
    with Pool() as workers:

        for r, pair_res in workers.imap(create_molecule_selection, contigs):
            for cell, sel in r.items():
                selection[cell].update(sel)
            pairings.update( pair_res )

    molecule = None


    with pysam.AlignmentFile(bam_path,threads=8) as alignments, \
        pysam.FastaFile(get_reference_path_from_bam(alignments)) as reference:
        with sorted_bam_file(output_path,alignments) as out:
            for molecule in MoleculeIterator(get_nearby_reads(),

                                             molecule_class=CHICMolecule,
                                             fragment_class = CHICFragment,
                                             molecule_class_args = {'reference':reference}
                                            ):

                molecule.write_pysam(out, consensus=True,consensus_read_callback=crc, consensus_read_callback_kwargs={'molecule':molecule})
