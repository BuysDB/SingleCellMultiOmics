#!/usr/bin/env python
# -*- coding: utf-8 -*-
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysamiterators
import pysam
import argparse
import singlecellmultiomics.bamProcessing.bamFunctions as bf
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
        '-head',
        type=int,
        help='Tabulate the first N valid molecules')
    argparser.add_argument('-minmq', type=int, default=50)
    argparser.add_argument(
        '-contig',
        type=str,
        help='contig to run on, all when not specified')
    argparser.add_argument('-method', type=str, default='nla')
    argparser.add_argument(
        '-fmt',
        type=str,
        default='table',
        help="output format (options are: 'bed' or 'table')")
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
    taps = singlecellmultiomics.molecule.TAPS(reference=reference)

    molecule_class_args = {
        'reference': reference,
        'taps': taps,
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

    try:
        for i, molecule in enumerate(singlecellmultiomics.molecule.MoleculeIterator(
            alignments=alignments,
            molecule_class=molecule_class,
            yield_invalid=(output is not None),
            fragment_class=fragment_class,
            fragment_class_args={'umi_hamming_distance': 1},

                molecule_class_args=molecule_class_args,
                contig=args.contig)):

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

            for (chromosome, location), call in molecule.methylation_call_dict.items():
                if call['context'] == '.':  # Only print calls concerning C's
                    continue

                # Skip non-selected contexts
                if contexts is not None and call['context'] not in contexts:
                    continue

                if args.fmt == "table":
                    print(
                        f"{molecule.sample}{args.moleculeNameSep}{i}{args.moleculeNameSep}{molecule.umi}{args.moleculeNameSep}{molecule.get_strand_repr()}\t{chromosome}\t{location+1}\t{call['context']}")
                elif args.fmt == "bed":
                    sample = molecule.sample.split("_")[-1]
                    print(
                        f'{chromosome}\t{location}\t{location+1}\t{sample}\t1\t{molecule.get_strand_repr()}')

    except (KeyboardInterrupt, BrokenPipeError) as e:
        pass

    if output is not None:
        finish_bam(output, args, temp_out)
