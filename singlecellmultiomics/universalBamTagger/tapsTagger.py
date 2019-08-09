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
import singlecellmultiomics.bamProcessing.bamFunctions as bf

class Fraction:
    def __init__(self):
        self.values = [0,0]
    def __setitem__(self, key, value):
        self.values[key] = value
    def __getitem__(self, key):
        return self.values[key]

    def __float__(self):
        if sum(self.values)==0:
            return np.nan
        return self.values[1]/(self.values[1]+self.values[0])

if __name__=='__main__':


    argparser = argparse.ArgumentParser(
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     description='Add methylation information to BAM file')

    argparser.add_argument('alignmentfile',  type=str)
    argparser.add_argument('-o',  type=str, help="output BAM path", required=True)
    argparser.add_argument('-table',  type=str, help="output table alias", required=False)

    argparser.add_argument('-bin_size',  type=int, default=250_000)
    argparser.add_argument('-sliding_increment',  type=int, default=50_000)


    argparser.add_argument('-ref',  type=str, required=False, help='path to reference fasta file ')
    argparser.add_argument('-min_mq',  default=20, type=int, help='Min mapping qual')
    argparser.add_argument('-uhd',  default=1, type=int, help='Umi hamming distance')
    argparser.add_argument('-head',  type=int)
    argparser.add_argument('-stranded',  action='store_true')
    argparser.add_argument('-contig',  type=str,help='contig to run on')

    argparser.add_argument('-samples',  type=str,help='Samples to select, separate with comma. For example CellA,CellC,CellZ', default=None)
    argparser.add_argument('-context',  type=str,help='Contexts to select, separate with comma. For example Z,H,X', default=None)
    args = argparser.parse_args()

    samples = None if args.samples is None else set(args.samples.split(','))
    contexts = None if args.context is None else set(
        [x.upper() for x in  args.context.split(',')] +
        [x.lower() for x in  args.context.split(',')])

    alignments= pysam.AlignmentFile(args.alignmentfile)
    # Auto detect reference:
    if args.ref is None:
        args.ref = bf.get_reference_from_pysam_alignmentFile(alignments)
        if args.ref is None:
            raise ValueError("Supply reference, -ref")

    if args.contig is None:
        # Create jobs for all chromosomes:
        temp_prefix = os.path.abspath( os.path.dirname(args.o) )+ '/' + str(uuid.uuid4())
        hold_merge=[]
        for chrom in alignments.references:
            if chrom.startswith('KN') or chrom.startswith('KZ') or chrom.startswith('chrUn') or chrom.endswith('_random') or 'ERCC' in chrom:
                continue
            temp_bam_path = f'{temp_prefix}_{chrom}.bam'
            arguments = " ".join([x for x in sys.argv if not x==args.o in x and x!='-o'])  + f" -contig {chrom} -o {temp_bam_path}"
            job = f'TAPS_{str(uuid.uuid4())}'
            os.system( f'submission.py' + f' -y --py36 -time 50 -t 1 -m 8 -N {job} " {arguments};"' )
            hold_merge.append(job)

        hold =  ','.join(hold_merge)
        os.system( f'submission.py' + f' -y --py36 -time 50 -t 1 -m 8 -N {job} -hold {hold} " samtools merge {args.o} {temp_prefix}*.bam; samtools index {args.o}; rm {temp_prefix}*.ba*"' )
        exit()

    reference = pysamiterators.iterators.CachedFasta( pysam.FastaFile(args.ref) )
    taps = singlecellmultiomics.molecule.TAPS(reference=reference)
    temp_out = f'{args.o}.temp.out.bam'

    # Obtain contig sizes:
    ref_lengths = {r:alignments.get_reference_length(r) for r in alignments.references}

    # Methylation dictionary: site->cell->value
    binned_data = collections.defaultdict(lambda: collections.defaultdict(Fraction))
    cell_count=collections.Counter()

    with pysam.AlignmentFile(temp_out , "wb",header=alignments.header) as output:
        for i,molecule in  enumerate(
            singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=alignments,
                    moleculeClass=singlecellmultiomics.molecule.TAPSNlaIIIMolecule,
                    fragmentClass=singlecellmultiomics.fragment.NLAIIIFragment,
                    fragment_class_args={'umi_hamming_distance':args.uhd},
                    yield_invalid=True,
                    molecule_class_args={
                       # 'features':transcriptome_features,
                        'reference':reference,
                        'site_has_to_be_mapped':True,
                        #radius in order to capture spliced transcripts
                        'taps':taps,
                        'min_max_mapping_quality':args.min_mq
                    },
                    contig=contig
            )):
            if i>args.head:
                break

            if samples is not None and molecule.sample not in samples:
                molecule.set_rejection_reason('sample_not_selected')
                if output is not None:
                    molecule.write_pysam(output)
                continue

            if not molecule.is_valid():
                molecule.set_meta('RF','rejected_molecule')
                molecule.write_tags()
                molecule.write_pysam(output)
                continue


            got_context_hit = False
            methylated_hits = 0
            unmethylated_hits = 0
            for (chromosome, location),call in molecule.methylation_call_dict.items():
                if call=='.': # Only use calls concerning C's
                    continue
                # Skip non-selected contexts
                if contexts is not None and call not in contexts:
                    continue
                got_context_hit+=1

                if call.isupper():
                    methylated_hits += 1
                else:
                    unmethylated_hits += 1
                if args.table is not None:
                    for binIdx in singlecellmultiomics.utils.coordinate_to_bins(location, args.bin_size, args.sliding_increment):
                        bin_start, bin_end = binIdx
                        if bin_start<0 or bin_end>ref_lengths[molecule.chromosome]:
                            continue

                        if args.stranded:
                            binned_data[(chromosome, molecule.get_strand_repr(), binIdx)][molecule.get_sample()][call.isupper()]+=1
                            cell_count[molecule.get_sample()]+=1
                        else:
                            binned_data[(chromosome, binIdx)][molecule.get_sample()][call.isupper()]+=1
                            cell_count[molecule.get_sample()]+=1

            molecule.set_meta('ME',methylated_hits)
            molecule.set_meta('um',unmethylated_hits)
            molecule.write_tags()
            molecule.set_meta('RF','accepted_molecule')
            molecule.write_pysam(output)

    if args.table is not None:
        # Write raw counts:
        df = pd.DataFrame(
            {loc:{sample: binned_data[loc][sample][0] for sample in binned_data[loc] } for loc in binned_data})
        df.to_pickle(f'{args.table}_unmethylated_{args.contig}.pickle.gz')
        df.to_csv(f'{args.table}_unmethylated_{args.contig}.csv')
        del df

        df = pd.DataFrame(
            {loc:{sample: binned_data[loc][sample][1] for sample in binned_data[loc] } for loc in binned_data})
        df.to_pickle(f'{args.table}_methylated_{args.contig}.pickle.gz')
        df.to_csv(f'{args.table}_methylated_{args.contig}.csv')
        del df

        #cast all fractions to float
        for loc in binned_data:
            for sample in binned_data[loc]:
                binned_data[loc][sample] = float(binned_data[loc][sample])
        df = pd.DataFrame(binned_data)
        del binned_data
        if args.contig:
            df.to_pickle(f'{args.table}_ratio_{args.contig}.pickle.gz')
            df.to_csv(f'{args.table}_ratio_{args.contig}.csv')
        else:
            df.to_pickle(f'{args.table}_ratio.pickle.gz')
            df.to_csv(f'{args.table}_ratio.csv')

    # Sort and index
    # Perform a reheading, sort and index
    cmd = f"""samtools sort {temp_out} > {args.o}; samtools index {args.o};
    rm {temp_out};
    """
    os.system(cmd)
