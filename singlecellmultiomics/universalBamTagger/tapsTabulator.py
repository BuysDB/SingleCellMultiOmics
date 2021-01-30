#!/usr/bin/env python
# -*- coding: utf-8 -*-
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysamiterators
import pysam
import argparse
import singlecellmultiomics.bamProcessing.bamFunctions as bf
from singlecellmultiomics.features import FeatureContainer
import os


def finish_bam(output, args, temp_out):
    output.close()
    # Sort and index
    # Perform a reheading, sort and index
    cmd = f"""samtools sort {temp_out} > {args.bamout}; samtools index {args.bamout};
    rm {temp_out};
    """
    os.system(cmd)


# sy zcat chr18 | sort -S 4G -k 3,3 -n -T ./TMP | gzip >
# chr18_sorted.gzip" -m 5

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract TAPS methylation calls from BAM file')

    argparser.add_argument('alignmentfile', type=str)
    argparser.add_argument(
        '-ref',
        type=str,
        help='path to reference fasta file, auto detected from bamfile')
    argparser.add_argument(
        '--every_fragment_as_molecule',
        action='store_true',
        help='Assign every fragment as a molecule, this effectively disables UMI deduplication and consensus calls based on multiple fragments')
    argparser.add_argument(
        '-head',
        type=int,
        help='Tabulate the first N valid molecules')


    argparser.add_argument(
        '-features',
        type=str,
        help='Annotate cut locations with these features, path to gtf file')


    argparser.add_argument(
        '-min_phred_score',
        type=int,
        help='Do not call methylation for bases with a phred score lower than "min_phred_score"',default=None)


    argparser.add_argument(
        '-dove_R1_distance',
        type=int,
        help='Do not call methylation N bases form the end of R1',default=8)

    argparser.add_argument(
        '-dove_R2_distance',
        type=int,
        help='Do not call methylation N bases form the end of R2',default=8)

    argparser.add_argument(
        '-skip_last_n_cycles_R1',
        type=int,
        help='Do not call methylation N bases form the end of R1',default=5)



    argparser.add_argument(
        '-skip_first_n_cycles_R1',
        type=int,
        help='Do not call methylation N bases form the start of R1',default=5)

    argparser.add_argument(
        '-skip_last_n_cycles_R2',
        type=int,
        help='Do not call methylation N bases form the end of R2',default=5)

    argparser.add_argument(
        '-skip_first_n_cycles_R2',
        type=int,
        help='Do not call methylation N bases form the start of R2',default=5)


    argparser.add_argument('-minmq', type=int, default=50)
    argparser.add_argument(
        '-contig',
        type=str,
        help='contig to run on, all when not specified')
    argparser.add_argument('-method', type=str, default='nla', help='nla, or chic')

    argparser.add_argument(
        '-fmt',
        type=str,
        default='table',
        help="""output format (options are: 'bed' or 'table' or 'table_more'), Columns are:
        table:
        sample, cut_site, molecule_span, umi, strand, chromosome, location (1 based), context
        long format (slower):
        sample, cut_site, molecule_span, mean_R1_cyle, mean_R2_cycle, mean_phred_score, umi, strand, chromosome, location (1 based), context
         """)
    argparser.add_argument(
        '-moleculeNameSep',
        type=str,
        help='Separator to use in molecule name',
        default=':')
    argparser.add_argument(
        '-samples',
        type=str,
        help='Samples to select, separate with comma. For example CellA,CellC,CellZ',
        default=None)
    argparser.add_argument(
        '-context',
        type=str,
        help='Contexts to select, separate with comma. For example Z,H,X',
        default=None)
    argparser.add_argument(
        '-bamout',
        type=str,
        help="optional (tagged) output BAM path")
    argparser.add_argument(
        '--allow_single_end',
        type=str,
        help='Allow single end reads')

    args = argparser.parse_args()
    alignments = pysam.AlignmentFile(args.alignmentfile)

    if args.bamout is not None:
        temp_out = f'{args.bamout}.unsorted.bam'
        output = pysam.AlignmentFile(temp_out, "wb", header=alignments.header)
    else:
        output = None

    samples = None if args.samples is None else set(args.samples.split(','))
    contexts = None if args.context is None else set(
        [x.upper() for x in args.context.split(',')] +
        [x.lower() for x in args.context.split(',')])

    if args.ref is None:
        args.ref = bf.get_reference_from_pysam_alignmentFile(alignments)
        if args.ref is None:
            raise ValueError("Supply reference, -ref")

    reference = pysamiterators.iterators.CachedFasta(pysam.FastaFile(args.ref))
    taps = singlecellmultiomics.molecule.TAPS()

    if args.features is not None:
        features = FeatureContainer()
        features.loadGTF(args.features,thirdOnly='gene',store_all=True)
    else:
        features = None

    molecule_class_args = {
        'reference': reference,
        'taps': taps,
        'taps_strand':'R',
        'min_max_mapping_quality': args.minmq
    }
    if args.method == 'nla':
        molecule_class = singlecellmultiomics.molecule.TAPSNlaIIIMolecule
        fragment_class = singlecellmultiomics.fragment.NlaIIIFragment
        molecule_class_args.update({'site_has_to_be_mapped': True})
    elif args.method == 'chic':
        molecule_class = singlecellmultiomics.molecule.TAPSCHICMolecule
        fragment_class = singlecellmultiomics.fragment.CHICFragment
    else:
        raise ValueError("Supply 'nla' or 'chic' for -method")


    molecule_class_args['methylation_consensus_kwargs'] = {

        'skip_first_n_cycles_R1':args.skip_first_n_cycles_R1,
        'skip_first_n_cycles_R2':args.skip_first_n_cycles_R2,
        'dove_R1_distance':args.dove_R1_distance,
        'dove_R2_distance':args.dove_R2_distance,
        'skip_last_n_cycles_R1':args.skip_last_n_cycles_R1,
        'skip_last_n_cycles_R2':args.skip_last_n_cycles_R2,
        'min_phred_score':args.min_phred_score
        }


    if args.allow_single_end:
        # Single end base calls are "unsafe", allow them :
        molecule_class_args['allow_unsafe_base_calls'] = True


    s = args.moleculeNameSep
    try:
        for i, molecule in enumerate(singlecellmultiomics.molecule.MoleculeIterator(
            alignments=alignments,
            molecule_class=molecule_class,
            yield_invalid=(output is not None),
            every_fragment_as_molecule=args.every_fragment_as_molecule,
            fragment_class=fragment_class,
            fragment_class_args={'umi_hamming_distance': 1,
                                 'no_umi_cigar_processing':False,

                                 },

                molecule_class_args=molecule_class_args,
                contig=args.contig)):

            molecule.set_meta('mi',i)
            if args.head and (i - 1) >= args.head:
                break

            if not molecule.is_valid(set_rejection_reasons=True):
                if output is not None:
                    molecule.write_pysam(output)
                continue

            # Skip sample if not selected
            if samples is not None and molecule.sample not in samples:
                molecule.set_rejection_reason('sample_not_selected')
                if output is not None:
                    molecule.write_pysam(output)
                continue

            if args.fmt == "table_more":
                consensus = molecule.get_consensus()

            cut_contig,CUT_SITE,cut_strand = molecule.get_cut_site()


            if features is not None:
                features.findFeaturesAt(cut_contig, CUT_SITE)
                genes = []
                for hit_start,hit_end,gene,strand,meta in features.findFeaturesAt(cut_contig, CUT_SITE):
                    g = dict(meta).get('gene_name',gene)
                    genes.append(g)
                additional=f"\t{','.join(genes)}"

            else:
                additional = ""


            for (chromosome, location), call in molecule.methylation_call_dict.items():
                if call['context'] == '.':  # Only print calls concerning C's
                    continue

                # Skip non-selected contexts
                if contexts is not None and call['context'] not in contexts:
                    continue

                if args.fmt == "table":

                    print(
                        f"{molecule.sample}{s}{i}{s}{CUT_SITE}{s}{molecule.estimated_max_length}{s}{molecule.umi}{s}{molecule.get_strand_repr()}\t{chromosome}\t{location+1}\t{call['context']}\t{molecule.ligation_motif}{additional}")

                elif args.fmt == "table_more":

                    base = consensus[chromosome, location]
                    bq = molecule.get_mean_base_quality(chromosome, location, base=base)
                    R1_cycle, R2_cycle = molecule.get_mean_cycle(chromosome, location, base=base)

                    print(
                        f"{molecule.sample}{s}{i}{s}{molecule.get_cut_site()[1]}{s}{molecule.span_len}{s}{R1_cycle:.0f}{s}{R2_cycle:.0f}{s}{bq:.0f}{s}{molecule.umi}{s}{molecule.get_strand_repr()}\t{chromosome}\t{location+1}\t{call['context']}")
                elif args.fmt == "bed":
                    sample = molecule.sample.split("_")[-1]
                    print(
                        f'{chromosome}\t{location}\t{location+1}\t{sample}\t1\t{molecule.get_strand_repr()}')

            if output is not None:
                molecule.write_tags()
                molecule.write_pysam(output)


    except (KeyboardInterrupt, BrokenPipeError) as e:
        pass

    if output is not None:
        finish_bam(output, args, temp_out)
