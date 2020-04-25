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
import singlecellmultiomics.features
import colorama
import numpy as np


class Fraction:
    def __init__(self):
        self.values = [0, 0]

    def __setitem__(self, key, value):
        self.values[key] = value

    def __getitem__(self, key):
        return self.values[key]

    def __float__(self):
        if sum(self.values) == 0:
            return np.nan
        return self.values[1] / (self.values[1] + self.values[0])


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Add methylation information to BAM file')

    argparser.add_argument('alignmentfile', type=str)
    argparser.add_argument(
        '-o',
        type=str,
        help="output BAM path",
        required=True)
    argparser.add_argument(
        '-bed',
        type=str,
        help="output bed file (base-level methylation status)",
        default=None,
        required=False)
    # these options need to go <- VB
    argparser.add_argument(
        '-table',
        type=str,
        help="output table alias",
        default=None,
        required=False)
    argparser.add_argument('-bin_size', type=int, default=250_000)
    argparser.add_argument('-sliding_increment', type=int, default=50_000)
    #
    argparser.add_argument(
        '-ref',
        type=str,
        required=False,
        help='path to reference fasta file ')
    argparser.add_argument(
        '-min_mq',
        default=20,
        type=int,
        help='Min mapping qual')
    argparser.add_argument('-uhd', default=1, type=int,
                           help='Umi hamming distance')
    argparser.add_argument(
        '-mem',
        default=40,
        type=int,
        help='Memory used per job')
    argparser.add_argument(
        '-time',
        default=52,
        type=int,
        help='Time requested per job')
    argparser.add_argument('-head', type=int)
    argparser.add_argument('-contig', type=str, help='contig to run on')
    argparser.add_argument('-method', type=str, help='nla or chic')
    argparser.add_argument(
        '-samples',
        type=str,
        help='Samples to select, separate with comma. For example CellA,CellC,CellZ',
        default=None)
    argparser.add_argument(
        '-context',
        type=str,
        help='Context to select, separate with comma. For example Z,H,X',
        default=None)

    argparser.add_argument('--stranded', action='store_true')
    argparser.add_argument(
        '--cluster',
        action='store_true',
        help='split by chromosomes and submit the job on cluster')
    argparser.add_argument(
        '--no_sort_index',
        action='store_true',
        help='do not sort and index the output bam')

    # Transcriptome splitting mode
    tr = argparser.add_argument_group('transcriptome specific settings')
    tr.add_argument(
        '--transcriptome',
        action='store_true',
        help='Label transcripts, requires exons and introns')
    tr.add_argument('-exons', type=str, help='Exon GTF file')
    tr.add_argument(
        '-introns',
        type=str,
        help='Intron GTF file, use exonGTF_to_intronGTF.py to create this file')
    tr.add_argument(
        '-recovery_umi_pool_radius',
        type=int,
        help='BP radius. When assigning transcripts without NLA site found, use this radius for molecule pooling',
        default=4)
    args = argparser.parse_args()

    samples = None if args.samples is None else set(args.samples.split(','))
    contexts = None if args.context is None else set(args.context.split(','))
#    set(
#        [x.upper() for x in  args.context.split(',')] +
#        [x.lower() for x in  args.context.split(',')])

    alignments = pysam.AlignmentFile(args.alignmentfile)
    # Auto detect reference:
    if args.ref is None:
        args.ref = bf.get_reference_from_pysam_alignmentFile(alignments)
        if args.ref is None:
            raise ValueError("Supply reference, -ref")

    if args.transcriptome:
        print(
            colorama.Style.BRIGHT +
            'Running in transcriptome recovery mode' +
            colorama.Style.RESET_ALL)
        if args.exons is None or args.introns is None:
            raise ValueError("Please supply both intron and exon GTF files")

    if args.cluster:
        if args.contig is None:
            # Create jobs for all chromosomes:
            temp_prefix = os.path.abspath(
                os.path.dirname(args.o)) + '/' + str(uuid.uuid4())
            hold_merge = []
            for chrom in alignments.references:
                if chrom.startswith('KN') or chrom.startswith('KZ') or chrom.startswith(
                        'chrUn') or chrom.endswith('_random') or 'ERCC' in chrom:
                    continue
                temp_bam_path = f'{temp_prefix}_{chrom}.bam'
                temp_bed_path = f'{temp_prefix}_{chrom}.bed'

                arguments = " ".join([x for x in sys.argv if not x == args.o and x != '-o']) + \
                    f" -contig {chrom} -o {temp_bam_path} -bed {temp_bed_path}"

                job = f'TAPS_{str(uuid.uuid4())}'
                os.system(
                    f'submission.py --silent' +
                    f' -y --py36 -time {args.time} -t 1 -m {args.mem} -N {job} " {arguments};"')
                hold_merge.append(job)

            hold = ','.join(hold_merge)
            os.system(
                f'submission.py --silent' +
                f' -y --py36 -time {args.time} -t 1 -m 10 -N {job} -hold {hold} " samtools merge {args.o} {temp_prefix}*.bam; samtools index {args.o}; rm {temp_prefix}*.ba*; cat {temp_prefix}*.bed > {args.bed}; rm {temp_prefix}*.bed"')
            exit()

    reference = pysamiterators.iterators.CachedFasta(pysam.FastaFile(args.ref))
    taps = singlecellmultiomics.molecule.TAPS(reference=reference)
    temp_out = f'{args.o}.temp.out.bam'

    # Obtain contig sizes:
    ref_lengths = {r: alignments.get_reference_length(
        r) for r in alignments.references}

    # Methylation dictionary: site->cell->value
    binned_data = collections.defaultdict(
        lambda: collections.defaultdict(Fraction))
    cell_count = collections.Counter()

    # Define molecule class arguments
    molecule_class_args = {
        'reference': reference,
        'taps': taps,
        'min_max_mapping_quality': args.min_mq
    }

    fragment_class_args = {'umi_hamming_distance': args.uhd}

    # transcriptome mode specific arguments: ####
    if args.transcriptome:

        transcriptome_features = singlecellmultiomics.features.FeatureContainer()
        print("Loading exons", end='\r')
        transcriptome_features.loadGTF(
            args.exons,
            select_feature_type=['exon'],
            identifierFields=(
                'exon_id',
                'gene_id'),
            store_all=True,
            contig=args.contig,
            head=None)

        print("Loading introns", end='\r')
        transcriptome_features.loadGTF(
            args.introns,
            select_feature_type=['intron'],
            identifierFields=['transcript_id'],
            store_all=True,
            contig=args.contig,
            head=None)
        print("All features loaded")

        rejected_reads = []  # Store all rejected, potential transcript reads

        # Add more molecule class arguments
        molecule_class_args.update({
            'features': transcriptome_features
        })

    # Method specific arguments
    if args.method == 'nla':
        molecule_class_args.update({'site_has_to_be_mapped': True})
    elif args.method == 'chic':
        fragment_class_args.update({'invert_strand': True})

    if args.transcriptome:
        if args.method == 'nla':
            molecule_class = singlecellmultiomics.molecule.AnnotatedTAPSNlaIIIMolecule
            fragment_class = singlecellmultiomics.fragment.NlaIIIFragment
        elif args.method == 'chic':
            molecule_class = singlecellmultiomics.molecule.AnnotatedTAPSCHICMolecule
            fragment_class = singlecellmultiomics.fragment.CHICFragment
        else:
            raise ValueError("Supply 'nla' or 'chic' for -method")
    else:
        if args.method == 'nla':
            molecule_class = singlecellmultiomics.molecule.TAPSNlaIIIMolecule
            fragment_class = singlecellmultiomics.fragment.NlaIIIFragment
        elif args.method == 'chic':
            molecule_class = singlecellmultiomics.molecule.TAPSCHICMolecule
            fragment_class = singlecellmultiomics.fragment.CHICFragment
        else:
            raise ValueError("Supply 'nla' or 'chic' for -method")

    ###############################################
    statistics = collections.defaultdict(collections.Counter)
    mcs = collections.Counter()  # methylation calls seen
    print(
        colorama.Style.BRIGHT +
        "Running TAPS tagging" +
        colorama.Style.RESET_ALL)
    if args.bed is not None:
        bed = open(args.bed, "w")
    with pysam.AlignmentFile(temp_out, "wb", header=alignments.header) as output:
        for i, molecule in enumerate(
                singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=alignments,
                    molecule_class=molecule_class,
                    fragment_class=fragment_class,
                    fragment_class_args=fragment_class_args,
                    yield_invalid=True,
                    molecule_class_args=molecule_class_args,
                    contig=args.contig
                )):

            if args.head is not None and i >= args.head:
                print(
                    colorama.Style.BRIGHT +
                    colorama.Fore.RED +
                    f"Head was supplied, stopped at {i} molecules" +
                    colorama.Style.RESET_ALL)
                break
            statistics['Input']['molecules'] += 1
            statistics['Input']['fragments'] += len(molecule)

            # Set (chromosome) unique identifier
            molecule.set_meta('mi', f'NLA_{i}')
            if args.transcriptome:
                molecule.set_intron_exon_features()

            if samples is not None and molecule.sample not in samples:
                molecule.set_rejection_reason('sample_not_selected')
                if output is not None:
                    molecule.write_pysam(output)
                continue

            if args.transcriptome:
                if not molecule.is_valid():
                    if molecule.is_multimapped() or molecule.get_max_mapping_qual() < args.min_mq:
                        molecule.set_meta('RF', 'rejected_molecule_mq')
                        molecule.write_tags()
                        molecule.write_pysam(output)
                        statistics['Filtering']['low mapping quality'] += 1
                        statistics['Filtering']['rejected'] += 1
                        continue

                    rejected_reads.append(molecule[0].reads)
                    continue
                statistics['Filtering'][f'valid {args.method} molecule'] += 1
                if len(molecule.junctions):
                    molecule.set_meta('RF', 'transcript_junction')
                    molecule.set_meta('dt', 'RNA')
                    statistics['Data type detection']['RNA because junction found'] += 1
                else:
                    if len(molecule.genes) == 0:
                        molecule.set_meta('dt', 'DNA')
                        statistics['Data type detection']['DNA not mapping to gene'] += 1
                    else:
                        # Check if NLA III sites are skipped...
                        if args.method == 'nla':
                            skipped = molecule.get_undigested_site_count()
                            if skipped == 0:
                                molecule.set_meta('dt', 'DNA')
                            else:
                                molecule.set_meta('dt', 'RNA or DNA')
                        else:
                            molecule.set_meta('dt', 'RNA or DNA')
            else:
                if not molecule.is_valid():
                    statistics['Filtering'][f'not valid {args.method}'] += 1
                    molecule.set_meta('RF', 'rejected_molecule')
                    molecule.write_tags()
                    molecule.write_pysam(output)
                    continue
                statistics['Filtering'][f'valid {args.method} molecule'] += 1
                molecule.set_meta('RF', 'accepted_molecule')

            got_context_hit = False
            methylated_hits = 0
            unmethylated_hits = 0
            readString = []
            genomeString = []
            for (chromosome, location), call in molecule.methylation_call_dict.items():
                if call['context'] == '.':  # Only use calls concerning C's
                    continue

                got_context_hit += 1
                mcs[call['context']] += 1
                if call['context'].isupper():
                    methylated_hits += 1
                    readString.append(call['consensus'])
                    genomeString.append(call['reference_base'])
                else:
                    unmethylated_hits += 1
# NEED TO REMOVE THIS CODE ENTIRELY!! <- VB
                if args.table is not None:
                    for binIdx in singlecellmultiomics.utils.coordinate_to_bins(
                            location, args.bin_size, args.sliding_increment):
                        bin_start, bin_end = binIdx
                        if bin_start < 0 or bin_end > ref_lengths[molecule.chromosome]:
                            continue

                        if args.stranded:
                            binned_data[(chromosome, molecule.get_strand_repr(
                            ), binIdx)][molecule.get_sample()][call['context'].isupper()] += 1
                            cell_count[molecule.get_sample()] += 1
                        else:
                            binned_data[(chromosome, binIdx)][molecule.get_sample(
                            )][call['context'].isupper()] += 1
                            cell_count[molecule.get_sample()] += 1
###
                if args.bed is not None:
                    # Skip non-selected contexts only for table
                    if contexts is not None and call['context'] not in contexts:
                        continue
                    else:
                        # name = cell barcode + context
                        name = ":".join(
                            [molecule.sample.split("_")[-1], call['context']])
                        bed.write(
                            f'{chromosome}\t{location}\t{location+1}\t{name}\t1\t{molecule.get_strand_repr()}\n')

            refbase = '' if not genomeString else max(
                set(genomeString), key=genomeString.count)
            readbase = '' if not readString else max(
                set(readString), key=readString.count)

            readConversionString = None
            genomeConversionString = None
            # OT
            readConversionString = None
            genomeConversionString = None
            if readbase == 'T' and refbase == 'C' and molecule.get_strand() == 1:  # '+'
                readConversionString = 'CT'
                genomeConversionString = 'CT'
            # OB
            elif readbase == 'A' and refbase == 'G' and molecule.get_strand() == 0:  # '-'
                readConversionString = 'CT'
                genomeConversionString = 'GA'
            # CTOT
            elif readbase == 'A' and refbase == 'C' and molecule.get_strand() == 1:  # '+'
                readConversionString = 'GA'
                genomeConversionString = 'CT'
            # CTOB
            elif readbase == 'A' and refbase == 'G' and molecule.get_strand() == 0:  # '-'
                readConversionString = 'GA'
                genomeConversionString = 'GA'

            if readConversionString is not None:
                molecule.set_meta('XR', readConversionString)
            if genomeConversionString is not None:
                molecule.set_meta('XG', genomeConversionString)

            molecule.set_meta('ME', methylated_hits)
            molecule.set_meta('um', unmethylated_hits)
            statistics['Methylation']['methylated Cs'] += methylated_hits
            statistics['Methylation']['unmethylated Cs'] += unmethylated_hits
            molecule.write_tags()
            molecule.write_pysam(output)
# close bed
    if args.bed is not None:
        bed.close()

    if args.transcriptome:
        print(
            colorama.Style.BRIGHT +
            f"Running transcriptome recovery on {len(rejected_reads)} reads")
        for i, molecule in enumerate(
            singlecellmultiomics.molecule.MoleculeIterator(
                # plug in the possible_transcripts as read source
                alignments=rejected_reads,
                # Drop the TAPS and NLAIII checks
                molecule_class=singlecellmultiomics.molecule.FeatureAnnotatedMolecule,
                # Plain fragment, no NLAIII
                fragment_class=singlecellmultiomics.fragment.Fragment,
                fragment_class_args={
                    'umi_hamming_distance': args.uhd,
                    # this is the amount of bases R1 can shift to be assigned to the same molecule
                    'assignment_radius': args.recovery_umi_pool_radius
                },

                yield_invalid=True,
                molecule_class_args={
                        'features': transcriptome_features,
                    'reference': reference,
                            'min_max_mapping_quality': 20
                }
            )):
            if not molecule.is_valid():
                statistics['Filtering']['rejected at transcriptome recovery step'] += 1
                statistics['Filtering']['rejected'] += 1
                molecule.set_meta('RF', 'rejected_recovery_invalid')
                molecule.write_tags()
                molecule.write_pysam(output)
                continue

            molecule.set_meta('mi', f'TRAN_{i}')

            # Add gene annotations:
            molecule.annotate(0)
            molecule.set_intron_exon_features()
            if len(molecule.genes) == 0:
                molecule.set_meta('RF', 'rejected_recovery_no_gene')
                statistics['Filtering']['rejected_recovery_no_gene'] += 1
                molecule.write_tags()
                molecule.write_pysam(output)
                continue
            if len(molecule.junctions):
                molecule.set_meta('RF', 'recovered_transcript_junction')
                statistics['Filtering']['recovered_transcript_junction'] += 1
                statistics['Data type detection'][
                    f'RNA because junction found and no {args.method} site mapped'] += 1
            else:
                molecule.set_meta('RF', 'recovered_transcript_gene')
                statistics['Data type detection'][
                    f'RNA because gene found and no {args.method} site mapped'] += 1
                statistics['Filtering']['recovered_transcript_gene'] += 1
            molecule.set_meta('dt', 'RNA')
            statistics['Data type detection']['RNA'] += 1
            molecule.write_tags()
            molecule.write_pysam(output)

    # Show statistics:
    print(
        '\n' +
        colorama.Style.BRIGHT +
        'Statistics' +
        colorama.Style.RESET_ALL)
    for statistic_class in [
        'Input',
        'Filtering',
        'Data type detection',
            'Methylation']:
        print(f'{colorama.Style.BRIGHT} {statistic_class} {colorama.Style.RESET_ALL}')
        for statistic, value in statistics[statistic_class].most_common():
            print(f'  {statistic}\t{value}')

    print(f'{colorama.Style.BRIGHT} Methylation calls {colorama.Style.RESET_ALL}')
    for call, description in zip('zZxXhH',
                                 ['unmethylated C in CpG context (CG)',
                                  'methylated C in CpG context (CG)',
                                  'unmethylated C in CHG context ( C[ACT]G )',
                                  'methylated C in CHG context   ( C[ACT]G )',
                                  'unmethylated C in CHH context ( C[ACT][ACT] )',
                                  'methylated C in CHH context ( C[ACT][ACT] )'
                                  ]):

        if call.isupper():
            print(
                f'  {colorama.Style.BRIGHT}{call}{colorama.Style.RESET_ALL}\t{mcs[call]}',
                end='\t')
        else:
            print(f'  {call}\t{mcs[call]}', end='\t')
        print(f'{colorama.Style.DIM}{description}{colorama.Style.RESET_ALL}')

    print('\n')

    if args.table is not None:
        print(
            colorama.Style.BRIGHT +
            'Writing raw unmethylated counts' +
            colorama.Style.RESET_ALL)
        # Write raw counts:
        df = pd.DataFrame(
            {loc: {sample: binned_data[loc][sample][0] for sample in binned_data[loc]} for loc in binned_data})
        df.to_pickle(f'{args.table}_unmethylated_{args.contig}.pickle.gz')
        df.to_csv(f'{args.table}_unmethylated_{args.contig}.csv')
        del df

        print(
            colorama.Style.BRIGHT +
            'Writing raw methylated counts' +
            colorama.Style.RESET_ALL)
        df = pd.DataFrame(
            {loc: {sample: binned_data[loc][sample][1] for sample in binned_data[loc]} for loc in binned_data})
        df.to_pickle(f'{args.table}_methylated_{args.contig}.pickle.gz')
        df.to_csv(f'{args.table}_methylated_{args.contig}.csv')
        del df

        print(
            colorama.Style.BRIGHT +
            'Writing ratio tables' +
            colorama.Style.RESET_ALL)
        # cast all fractions to float
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

    print(
        colorama.Style.BRIGHT +
        'Sorting and indexing final file' +
        colorama.Style.RESET_ALL)
    # Sort and index
    # Perform a reheading, sort and index
    if not args.no_sort_index:
        cmd = f"""samtools sort {temp_out} > {args.o}; samtools index {args.o};
        rm {temp_out};
        """
    else:
        cmd = f"mv {temp_out} {args.o}"
    os.system(cmd)
    print("All done.")
