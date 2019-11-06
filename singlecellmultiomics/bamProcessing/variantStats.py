#!/usr/bin/env python3
import pysam
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import singlecellmultiomics.pyutils
import pysamiterators
import collections
import glob
import pickle
import pandas as pd
from colorama import Fore, Style
from singlecellmultiomics.bamProcessing.bamFunctions import sorted_bam_file, sort_and_index, get_reference_from_pysam_alignmentFile, add_readgroups_to_header, write_program_tag
from singlecellmultiomics.pyutils import meanOfCounter
import argparse

f'Please use an environment with python 3.6 or higher!'

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""
Known variant locations extraction tool
""")
argparser.add_argument('bamfiles', nargs='+')
argparser.add_argument('-ssnv', help="sSNV bed file", type=str, required=True)
argparser.add_argument('-gsnv', help="gSNV bed file", type=str, required=True)
argparser.add_argument('-reference', help="reference fasta file", type=str, required=True)
argparser.add_argument('-head', type=int)
args=  argparser.parse_args()
WINDOW_RADIUS = 250

paths = args.bamfiles

def obtain_variant_statistics(
    alignment_file_paths,
    cell_obs, statistics,cell_call_data,
    reference,
    chromosome,
    ssnv_position,gsnv_position,haplotype_scores,
    WINDOW_RADIUS, out

    ):


    sSNV_ref_base = reference.fetch(chromosome,ssnv_position,ssnv_position+1)
    gSNV_ref_base = reference.fetch(chromosome,gsnv_position,gsnv_position+1)

    window_molecules = []

    for pathi,path in enumerate(alignment_file_paths):
        with pysam.AlignmentFile(path) as alignments:

            for molecule_id,molecule in enumerate(
                singlecellmultiomics.molecule.MoleculeIterator(alignments,
                    fragment_class_args={
                        'umi_hamming_distance':1,

                    },
                    molecule_class_args={
                        'reference':reference
                    },
                   moleculeClass=singlecellmultiomics.molecule.NlaIIIMolecule,
                   fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
                    start = ssnv_position-WINDOW_RADIUS,
                    end = ssnv_position+ WINDOW_RADIUS,
                   contig=chromosome
                )):

                # For every molecule obtain the consensus from which to extract the gSNV and sSNV:
                try:
                    consensus = molecule.get_consensus()
                except Exception as e:
                    if str(e) == 'Could not extract any safe data from molecule':
                        statistics[(chromosome, ssnv_position)]['R2_unmapped'][True]+=1
                    continue

                # Extract the gSNV and sSNV:
                ssnv_state = consensus.get((chromosome, ssnv_position))
                gsnv_state = consensus.get((chromosome, gsnv_position) )
                window_molecules.append((molecule,ssnv_state,gsnv_state))

                # If both are none there is no information we can use.
                if ssnv_state is None and gsnv_state is None:
                    continue

                # Store the observation
                cell_obs[(chromosome, ssnv_position)][molecule.get_sample()][(ssnv_state,gsnv_state)] += 1

                # Store statistics
                statistics[(chromosome, ssnv_position)]['max_mapping_quality'][ molecule.get_max_mapping_qual() ] +=1
                statistics[(chromosome, ssnv_position)]['fragment_size'][molecule.get_safely_aligned_length()]+=1
                statistics[(chromosome, ssnv_position)]['ivt_dups'][ len(  molecule.get_rt_reactions() )] +=1
                statistics[(chromosome, ssnv_position)]['undigested'][ molecule.get_undigested_site_count()]+=1
                statistics[(chromosome, ssnv_position)]['reads'][ len(molecule) ] +=1
                statistics[(chromosome, ssnv_position)]['molecules'][1] +=1

                for operation, per_bp in molecule.get_alignment_stats().items():
                    statistics[(chromosome, ssnv_position)][operation][per_bp] +=1

                try:
                    statistics[(chromosome, ssnv_position)]['ssnv_ref_phred'][ molecule.get_mean_base_quality(chromosome, ssnv_position, sSNV_ref_base ) ] +=1
                except:
                    pass

                try:
                    statistics[(chromosome, ssnv_position)]['ssnv_alt_phred'][ molecule.get_mean_base_quality(chromosome, ssnv_position, not_base = sSNV_ref_base )]+=1
                except:
                    pass

                try:
                    statistics[(chromosome, ssnv_position)]['gsnv_ref_phred'][ molecule.get_mean_base_quality(chromosome, gsnv_position, gSNV_ref_base)]+=1
                except:
                    pass

                try:
                    statistics[(chromosome, ssnv_position)]['gsnv_alt_phred'][ molecule.get_mean_base_quality(chromosome, gsnv_position, not_base= gSNV_ref_base)]+=1
                except:
                    pass

    # After finishing iteration over all molecules assign haplotypes
    chrom, pos = chromosome, ssnv_position #
    obs_for_cells = cell_obs[(chrom, pos)]

    sSNV_alt_base = None
    gSNV_alt_base = None

    genotype_obs = collections.Counter()
    complete_genotype_obs = collections.Counter()
    sSNV_obs_phased = collections.Counter()
    gSNV_obs_phased = collections.Counter()

    sSNV_obs = collections.Counter()
    gSNV_obs = collections.Counter()

    for cell, cell_data in obs_for_cells.items():
        for ssnv,gsnv in cell_data:
            genotype_obs[(ssnv,gsnv)]+=1

            gSNV_obs[gsnv] +=1 # Only count these when the germline variant is detected
            sSNV_obs[ssnv] += 1

            if ssnv is not None and gsnv is not None:
                complete_genotype_obs[(ssnv,gsnv)]+=1
                gSNV_obs_phased[gsnv] +=1 # Only count these when the germline variant is detected
                sSNV_obs_phased[ssnv] += 1

    print(Style.BRIGHT + f'Genotype observations for variant {chrom}:{pos}' + Style.RESET_ALL)
    print('som\tgerm\tobs')
    for (ssnv,gsnv), obs  in complete_genotype_obs.most_common():
        print(f' {ssnv}\t{gsnv}\t{obs}')

    if len(complete_genotype_obs)<=2:
        print(f'not enough genotype observations for a variant call (<=2)')

    if len(complete_genotype_obs)>2:
        # Determine the sSNV alt:
        sSNV_alt_base=None
        for basecall, obs in sSNV_obs_phased.most_common():
            if basecall!=sSNV_ref_base:
                sSNV_alt_base = basecall
                break

        gSNV_alt_base=None
        for basecall, obs in gSNV_obs_phased.most_common():
            if basecall!=gSNV_ref_base:
                gSNV_alt_base = basecall
                break

        if sSNV_alt_base is None or gSNV_alt_base is None:
            # No phased alt base found ...
            return

        # Determine the phase (most common genotypes)
        print(f'sSNV alt is {sSNV_alt_base}')
        sSNV_phase = None
        wt_allele_gSNV= None
        for (sSNV_state, gSNV_state), obs in complete_genotype_obs.most_common():
            if sSNV_state == sSNV_alt_base:
                sSNV_phase = (sSNV_state, gSNV_state)
                phased_gSNV = gSNV_state
                if gSNV_state==gSNV_ref_base:
                    # the reference allele is alt
                    wt_allele_gSNV = gSNV_alt_base
                    snv_allele_gSNV = gSNV_ref_base
                else:
                    wt_allele_gSNV = gSNV_ref_base
                    snv_allele_gSNV = gSNV_alt_base
                    # the reference allele is ref

                break

        # The valid tuples are:
        valid_tuples = [ sSNV_phase,  # mutated
                        (sSNV_ref_base, snv_allele_gSNV), #wt
                        (sSNV_ref_base, wt_allele_gSNV)]

        # Count how often we found valid and invalid haplotypes
        valid = 0
        invalid = 0
        valid_var = 0
        invalid_var = 0

        for (ssnv,gsnv), tuple_obs in complete_genotype_obs.most_common():
            if ssnv == sSNV_alt_base: # variant:
                if (ssnv,gsnv) in valid_tuples:
                    valid_var+=tuple_obs
                else:
                    invalid_var+=tuple_obs

            if (ssnv,gsnv) in valid_tuples:
                valid += tuple_obs
            else:
                invalid +=tuple_obs

        phase_ratio = 0
        if valid_var+invalid_var > 0  :
            phase_ratio = valid_var/(valid_var+invalid_var)

        # Score Tuples with evidence for variant
        haplotype_scores[(chrom,pos)] = {
            'valid_tuples':valid,
            'invalid_tuples':invalid,
            'valid_var_tuples':valid_var,
            'invalid_var_tuples':invalid_var,
            'phasing_ratio' : phase_ratio,
            'gSNV_allelic_bias':gSNV_obs[gSNV_ref_base]/(gSNV_obs[gSNV_ref_base]+gSNV_obs[gSNV_alt_base])
        }


        print(f'Germline variant obs: {gSNV_ref_base} {gSNV_alt_base}')
        print(f'sSNV obs: {sSNV_ref_base} {sSNV_alt_base}')
        if sSNV_phase is not None:
            print(f'sSNV variant is phased with {phased_gSNV}')
        print(Style.BRIGHT +'Valid tuples:'+Style.RESET_ALL)
        for g,s in valid_tuples:
            print(f' {g}\t{s}')

        print(Style.BRIGHT +'Scores:'+Style.RESET_ALL)
        for name, obs in haplotype_scores[(chrom,pos)].items():
            print(f' {name}\t{obs}')

        # Create the cell call dictionary
        for cell, cell_data in obs_for_cells.items():
            for ssnv,gsnv in cell_data:
                if ssnv is None:
                    continue

                if gsnv is not None and not (ssnv,gsnv) in valid_tuples:
                    cell_call_data[(chrom, pos)][cell] =  -1 # invalid
                    continue

                if ssnv == sSNV_alt_base:
                    cell_call_data[(chrom, pos)][cell] = 1 # alternative
                    continue

                if ssnv == sSNV_ref_base and gsnv==phased_gSNV:
                    cell_call_data[(chrom, pos)][cell] = 0 # reference
                    continue

        # Annotate every molecule...

        for molecule_id,(m,ssnv_state, gsnv_state) in enumerate(window_molecules):
            m.set_meta('mi',molecule_id)
            if gsnv_state is None:
                m.set_meta('gv','?')
            else:
                m.set_meta('gv',gsnv_state)

            if ssnv_state is None:
                m.set_meta('sv','?')
            else:
                m.set_meta('sv',ssnv_state)


            if ssnv_state is None:
                m.set_meta('VD','NO_SNV_OVERLAP')
                continue

            if gsnv_state is not None and not (ssnv_state,gsnv_state) in valid_tuples:
                m.set_meta('VD','INVALID_PHASE')
                continue
            if ssnv_state == sSNV_alt_base:
                m.set_meta('VD','SNV_ALT')
                continue

            if ssnv_state == sSNV_ref_base and gsnv_state==phased_gSNV:
                m.set_meta('VD','SNV_REF')
                continue
            if gsnv_state!=phased_gSNV:
                m.set_meta('VD','UNINFORMATIVE_ALLELE')
                continue

            m.set_meta('VD','REJECTED')
        # write

        for m,ssnv_state, gsnv_state in window_molecules:
            m.write_tags()
            m.write_pysam(out)


# Load probed variants
probed_variants = {}
with open(args.ssnv) as s, \
     open(args.gsnv) as g:
    for i,(ssnv_line,gsnv_line) in enumerate(zip(s, g)):
        if ssnv_line.startswith('track name'):
            continue
        chrom, snv_pos, _ = ssnv_line.strip().split()
        _, gsnv_pos,__ = gsnv_line.strip().split()
        snv_pos, gsnv_pos = int(snv_pos), int(gsnv_pos)
        probed_variants[ (chrom,snv_pos) ]  = gsnv_pos

reference = pysamiterators.CachedFasta(pysam.FastaFile(args.reference))


cell_obs = collections.defaultdict(lambda: collections.defaultdict( collections.Counter) )
statistics = collections.defaultdict(lambda: collections.defaultdict( collections.Counter) )
cell_call_data =collections.defaultdict(dict) #location->cell->haplotype
haplotype_scores = {}


with sorted_bam_file('evidence.bam', origin_bam=pysam.AlignmentFile(paths[0]) ) as out:
    for variant_index, ((chromosome, ssnv_position),potential_gsnv_position) in enumerate(probed_variants.items()):

        obtain_variant_statistics(
            alignment_file_paths = paths,
            cell_obs = cell_obs,
            cell_call_data=cell_call_data,
            statistics = statistics,
            reference=reference,
            chromosome = chromosome,
            ssnv_position = ssnv_position,
            gsnv_position = potential_gsnv_position,
            WINDOW_RADIUS=WINDOW_RADIUS,
            haplotype_scores=haplotype_scores,
            out=out
        )

        if args.head and (variant_index>args.head-1):
            print(f'Stopping at variant {variant_index+1} because head was supplied ')
            break


lambda_free_dict = {}
for key, stats in statistics.items():
    lambda_free_dict[key] = {

        'mean_clip_pbp' : meanOfCounter(stats['clips_per_bp']),
        'mean_ins_pbp' : meanOfCounter(stats['inserts_per_bp']),
        'mean_del_pbp' : meanOfCounter(stats['deletions_per_bp']),
        'mean_matches_pbp' : meanOfCounter(stats['matches_per_bp']),
        'mean_alt_mapping_per_read' : meanOfCounter(stats['alt_per_read']),

        'ssnv_ref_phred' : meanOfCounter(stats['ssnv_ref_phred']),
        'ssnv_alt_phred' : meanOfCounter(stats['ssnv_alt_phred']),
        'gsnv_ref_phred' : meanOfCounter(stats['gsnv_ref_phred']),
        'gsnv_alt_phred' : meanOfCounter(stats['gsnv_alt_phred']),

        'mean_max_mapping_quality' : meanOfCounter(stats['max_mapping_quality']),
        'mean_ivt_dups' : meanOfCounter(stats['ivt_dups']),
        'mean_undigested' : meanOfCounter(stats['undigested']),
        'R2_unmapped' : stats['R2_unmapped'][True],


        'mean_fragment_size' : meanOfCounter(stats['fragment_size']),
        'mean_reads' : meanOfCounter(stats['reads']),
        'total_reads' : sum( (amount*frequency for amount, frequency in stats['reads'].most_common() )),
        'total_molecules' : sum( (amount*frequency for amount, frequency in stats['molecules'].most_common() ))

    }



print('Writing final site table')
site_stats = pd.DataFrame( lambda_free_dict ).T.join(pd.DataFrame(haplotype_scores).T)
site_stats.to_pickle('site_stats.pickle.gz')
site_stats.to_csv('site_stats.csv')

print('Writing final cell table')
cell_call_df = pd.DataFrame(cell_call_data)
cell_call_df.to_pickle('cell_calls.pickle.gz')
cell_call_df.to_csv('cell_calls.csv')
