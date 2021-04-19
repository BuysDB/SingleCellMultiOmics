#!/usr/bin/env python
# -*- coding: utf-8 -*-

from singlecellmultiomics.molecule import MoleculeIterator, CHICMolecule
from singlecellmultiomics.fragment import CHICFragment
from singlecellmultiomics.bamProcessing.bamFunctions import get_reference_path_from_bam, get_contigs
from collections import defaultdict
import pyBigWig
import pysam
from singlecellmultiomics.bamProcessing.bamBinCounts import get_binned_counts
import pandas as pd
import argparse
from singlecellmultiomics.bamProcessing import get_contig_sizes
from singlecellmultiomics.utils import is_autosome, pool_wrapper
import numpy as np
from multiprocessing import Pool
from more_itertools import windowed
from glob import glob
import pyBigWig
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import argparse
from itertools import product
import asyncio

def get_aligned_chic_len(molecule):

    for f in molecule:
        f.R2_primer_length = 6

    mlen = molecule.estimated_max_length
    if mlen is None:
        return None

    return mlen + 2 # 2 extra bases because 1 base is potentially lost on both sides

def _generate_molecule_coordinate_dict(args, max_fragment_size=800):

    bam_path,  sel_contig= args
    molecule_coordinate_dict = defaultdict(list) # cell->[ (start,end), (start,end) ..]
    with pysam.AlignmentFile(bam_path) as alignments:
       # pysam.FastaFile(get_reference_path_from_bam(alignments)) as reference:

        for molecule in MoleculeIterator(alignments,
                                         contig=sel_contig,
                                         molecule_class=CHICMolecule,
                                         fragment_class = CHICFragment,
                                         #molecule_class_args = {'reference':reference}
                                        ):



            contig, start, strand = molecule.get_cut_site()
            if not molecule.is_completely_matching:
                continue

            if molecule.get_mean_mapping_qual()<55:
                continue
            slen = get_aligned_chic_len(molecule)

            if slen is None or slen>max_fragment_size:
                continue

            if strand:
                end = start - slen
            else:
                end = start + slen


            molecule_coordinate_dict[molecule.sample].append( (start,end) )
    for sample in molecule_coordinate_dict:
        molecule_coordinate_dict[sample] =  sorted(molecule_coordinate_dict[sample])
    return sel_contig, molecule_coordinate_dict


def generate_molecule_coordinate_dict(bams, n_threads=None):

    molecule_coordinate_dicts = {}
    contigs = get_contigs(bams[0])

    with Pool(n_threads) as workers:
        for contig,d in workers.imap( _generate_molecule_coordinate_dict,
                                     product(bams,filter(is_autosome, contigs))):
            if not contig in molecule_coordinate_dicts:

                molecule_coordinate_dicts[contig] = d
            else:
                molecule_coordinate_dicts[contig].update( d )
            print(f'Finished {contig}')
    return molecule_coordinate_dicts


    #bin_size = 100_000

#distances_per_bin = defaultdict(lambda: defaultdict(list)  )

def contig_coordinate_dict_to_votes(d, n, min_total_cuts_per_cell=3, p_nuc_bin_size=5, max_linker_length=90, linker_vote_radius=25):
    #molecule_coordinate_dicts[contig] = d
    # n = size of contig
    vote_nucleosome = np.zeros(int( n/p_nuc_bin_size) )
    vote_linker = np.zeros(int( n/p_nuc_bin_size) )

    for cell in d.keys():
        if len(d[cell])<min_total_cuts_per_cell:
            continue
        for (start_a, end_a),(start_b,end_b) in windowed( d[cell],2):

            o_start_a = start_a

            #Check if the fragments are not overlapping:
            (start_a, end_a) = min(start_a, end_a), max(start_a, end_a)
            (start_b, end_b) = min(start_b, end_b), max(start_b, end_b)

            # The locations of the nucleosomes are constrained:
            # \\ start_a  -----  end_a \\     \\ start_b  ---- end_b \\
            if end_a - start_a >= 147 and end_a - start_a < (147*2+10): # it could contain a single nucleosome, but not more
                s = int(((start_a + 70 ))/p_nuc_bin_size)
                e = int(((end_a - 70))/p_nuc_bin_size)
                vote_nucleosome[ s:(e+1) ] += 1/(e-s)

            # Any starting point of a molecule is _always_ part of a linker
            # lets say += 25bp

            c=o_start_a
            s = int(((c-linker_vote_radius))/p_nuc_bin_size)
            e = int(((c+linker_vote_radius))/p_nuc_bin_size)+1
            vote_linker[s:(e+1)] += 1/(e-s)

            if end_a > start_b:
                continue

            if start_b - end_a  > max_linker_length: # The distance is larger than 90bp, which is a very long linker. Skip the linker vote
                continue

            s = int(((end_a))/p_nuc_bin_size)
            e = int(((start_b))/p_nuc_bin_size)


            vote_linker[s:(e+1)] += 1/(e-s+1)
    return vote_nucleosome, vote_linker


async def write_to_bw(handle, starts, ends,  values, contig,size=None):

    if size is not None:
        print( ends[ends>=size] )

    handle.addEntries(
        [contig]*len(starts), #Contig
        starts.astype(np.int64) , #Start
        ends= ends.astype(np.int64) ,  #end
        values=  values.astype(np.float64)
    )

async def coordinate_dicts_to_nucleosome_position_files(bams, molecule_coordinate_dicts, p_nuc_bin_size = 5, alias='npos', n_threads=None ):

    contigs = get_contigs(bams[0])
    sizes = get_contig_sizes(bams[0])


    # Create output bigwig file handles and worker pool
    linker_vote_write_path= f'{alias}_linkers.bw'
    nuc_vote_write_path= f'{alias}_nucleosomes.bw'
    merged_vote_write_path= f'{alias}_nucleosomes_min_linkers.bw'
    centroids = f'{alias}_nucleosome_centers.bed'


    with Pool(n_threads) as workers, \
        pyBigWig.open(linker_vote_write_path,'w') as linkers_out, \
        pyBigWig.open(nuc_vote_write_path,'w') as nuc_out, \
        pyBigWig.open(merged_vote_write_path,'w') as merged_out, \
        open(centroids,'w') as centroids_out:

        # Add headers for all output bigwigfiles
        for h in (linkers_out, nuc_out, merged_out):
            # Write header
            h.addHeader(list(zip(contigs, [sizes[c] for c in contigs])))


        # Obtain nucleosome and linker votes for each contig
        for ret_contig, (vote_nucleosome, vote_linker) in zip(
                    molecule_coordinate_dicts.keys(),
                    workers.imap(pool_wrapper,
            (
            (contig_coordinate_dict_to_votes,
            {
                'd':d,
                'n':sizes[contig],
                'min_total_cuts_per_cell':3,
                'p_nuc_bin_size':5,
                'max_linker_length':90,
                'linker_vote_radius':25
            })
                for contig, d in molecule_coordinate_dicts.items()
            ))):

            print(f'Writing votes for {ret_contig}')
            contig_len = sizes[ret_contig]

            smooth_linkers = gaussian_filter1d(vote_linker,2)
            smooth_nuc = gaussian_filter1d(vote_nucleosome,2)

            # Write nucleosome vote track:
            starts = np.array( np.linspace(0, contig_len-p_nuc_bin_size-1, len(vote_nucleosome)) )

            size = sizes[ret_contig]
            # Check if all numbers are preceding...
            d = np.diff(starts)
            d = d[d<1]
            if len(d):
                print(d)
            else:
                print(f'{len(starts)} Coordinates are in correct order')


            await asyncio.gather(
                write_to_bw(nuc_out, starts, starts+p_nuc_bin_size,  np.nan_to_num(smooth_nuc), ret_contig, size),
                write_to_bw(linkers_out, starts, starts+p_nuc_bin_size,  np.nan_to_num(smooth_linkers), ret_contig, size),
                write_to_bw(merged_out, starts, starts+p_nuc_bin_size,  np.nan_to_num(smooth_nuc - smooth_linkers), ret_contig, size=size)
            )


            # Use the find peak function to select properly spaced nucleosome positions.
            print(f'Writing centroids for {ret_contig}')
            mindist = 147
            min_dist_bins = int( mindist / p_nuc_bin_size )
            estimated_nucleosome_positions = find_peaks( gaussian_filter1d(smooth_nuc - smooth_linkers,2),distance=min_dist_bins)[0]

            for nucpos, score in zip( estimated_nucleosome_positions*p_nuc_bin_size, smooth_nuc[estimated_nucleosome_positions]):
                centroids_out.write(f'{ret_contig}\t{nucpos}\t{nucpos+p_nuc_bin_size}\t{score}\n')



    #contig_coordinate_dict_to_votes(d, n, min_total_cuts_per_cell=3, p_nuc_bin_size=5, max_linker_length=90, linker_vote_radius=25)


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Estimate nucleosome positions from scCHiC seq bam file(s)')

    argparser.add_argument('alignmentfiles', type=str, nargs='+')
    argparser.add_argument('-o', type=str, required=True, help='output prefix')

    argparser.add_argument('-bin_size', type=int,  default=3, help='Nucleosome position precision (bp), increasing this value increases memory consumption linearly')
    argparser.add_argument('-n_threads', type=int, help='Amount of worker threads')

    args = argparser.parse_args()

    async def main(args):
        # Run the analysis:
        print('Creating molecule coordinate database')
        d = generate_molecule_coordinate_dict(args.alignmentfiles, args.n_threads)
        print('Estimating nucleosome positions')
        await coordinate_dicts_to_nucleosome_position_files(args.alignmentfiles,
                    d,
                    p_nuc_bin_size = args.bin_size,
                    n_threads=args.n_threads,
                    alias=args.o)


    loop = asyncio.get_event_loop()
    loop.run_until_complete(main(args))

    #asyncio.run(main(args)) # python >= 3.7
