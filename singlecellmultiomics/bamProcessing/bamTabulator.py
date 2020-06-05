#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import pandas as pd
import singlecellmultiomics
import singlecellmultiomics.modularDemultiplexer
TagDefinitions = singlecellmultiomics.modularDemultiplexer.TagDefinitions
if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Tabulate a bam file to a file where every line corresponds to the features of a single read')
    argparser.add_argument(
        '-o',
        type=str,
        help="output csv path",
        required=False)
    argparser.add_argument('alignmentfile', type=str, nargs='?')
    argparser.add_argument('featureTags', type=str, nargs='?')
    argparser.add_argument(
        '-head',
        type=int,
        help='Run the algorithm only on the first N reads to check if the result looks like what you expect.')
    argparser.add_argument(
        '-e',
        action='store_true',
        help='Skip reads where any attribute is missing')
    argparser.add_argument(
        '--dedup',
        action='store_true',
        help='Count only the first occurence of a molecule. Requires RC tag or .duplicate bit to be set.')

    argparser.add_argument(
        '--noqcfail',
        action='store_true',
        help='Do not show qcfailed reads')

    argparser.add_argument(
        '--showtags',
        action='store_true',
        help='Show a list of commonly used tags, and tags present in your bam file')
    args = argparser.parse_args()

    if args.featureTags is None:
        args.o = None
        args.featureTags = None
        args.showtags = True

    if args.showtags or (args.alignmentfile is None and args.featureTags is None):
        # Find which tags are available in the file:
        head = 1000
        tagObs = collections.Counter()


        import colorama

        print(f'{colorama.Style.BRIGHT}Available attributes:{colorama.Style.RESET_ALL}')


        # Description from PYSAM: at https://pysam.readthedocs.io/en/latest/api.html
        attributes = [
            ('is_duplicate', 'true if optical or PCR duplicate' ),
            ('is_paired', 'true if read is paired in sequencing' ),

            ('is_proper_pair', 'true if read is mapped in a proper pair (reads are mapped within expected insert size)' ),
            ('is_qcfail', 'true if QC failure, rejected by tagger' ),

            ('is_read1', 'true if this is read1' ),
            ('is_read2', 'true if this is read2' ),
            ('is_reverse', 'true if read is mapped to reverse strand' ),
            ('is_secondary', 'true if not primary alignment' ),
            ('is_supplementary', 'true if this is a supplementary alignment' ),
            ('is_unmapped', 'true if read itself is unmapped' ),
            ('template_length', '' ),
            ('mapping_quality', 'mapping_quality' ),
            ('mate_is_reverse', 'true is read is mapped to reverse strand' ),
            ('mate_is_unmapped', '    true if the mate is unmapped' ),
            ('next_reference_id', 'Reference id of next mate' ),
            ('next_reference_name', 'Reference name (chromosome) of mate' ),
            ('query_alignment_end', 'end index of the aligned query portion of the sequence (0-based, exclusive)' ),
            ('query_alignment_length', 'length of the aligned query sequence.' ),

            ('template_length', 'length of the aligned query sequence.' ),
            ('query_alignment_qualities', """    aligned query sequence quality values (None if not present). These are the quality values that correspond to query, that is, they exclude qualities of soft clipped bases. This is equal to qual[qstart:qend]. Quality scores are returned as a python array of unsigned chars. Note that this is not the ASCII-encoded value typically seen in FASTQ or SAM formatted files. Thus, no offset of 33 needs to be subtracted.""" ),
            ('query_alignment_sequence', 'aligned portion of the read.' ),
            ('query_alignment_start', '' ),
            ('cigarstring', 'Alignment description' ),
            ('flag', 'properties flag' ),

        ]
        for attr,val in attributes:
            print(f'{colorama.Style.BRIGHT}{attr}{colorama.Style.RESET_ALL}\t{val}{colorama.Style.DIM}{colorama.Style.RESET_ALL}')

        if args.alignmentfile is not None:

            with pysam.AlignmentFile(args.alignmentfile, ignore_truncation=True) as f:
                for i, read in enumerate(f):
                    tagObs += collections.Counter([k for k,
                                                   v in read.get_tags(with_value_type=False)])
                    if i == (head - 1):
                        break

            print(
                f'{colorama.Style.BRIGHT}Tags seen in the supplied bam file(s):{colorama.Style.RESET_ALL}')
            for observedTag in tagObs:
                tag = observedTag
                if observedTag in TagDefinitions:
                    t = TagDefinitions[observedTag]
                    humanName = t.humanName
                    isPhred = t.isPhred
                else:
                    t = observedTag
                    isPhred = False
                    humanName = f'{colorama.Style.RESET_ALL}<No information available>'

                allReadsHaveTag = (tagObs[tag] == head)
                color = colorama.Fore.GREEN if allReadsHaveTag else colorama.Fore.YELLOW
                print(
                    f'{color}{ colorama.Style.BRIGHT}{tag}{colorama.Style.RESET_ALL}\t{color+humanName}\t{colorama.Style.DIM}{"PHRED" if isPhred else ""}{colorama.Style.RESET_ALL}' +
                    f'\t{"" if allReadsHaveTag else "Not all reads have this tag"}')

        print(f'{colorama.Style.BRIGHT}\nSCMO tag list:{colorama.Style.RESET_ALL}')
        for tag, t in TagDefinitions.items():
            print(f'{colorama.Style.BRIGHT}{tag}{colorama.Style.RESET_ALL}\t{t.humanName}\t{colorama.Style.DIM}{"PHRED" if t.isPhred else ""}{colorama.Style.RESET_ALL}')
        exit()

    if args.alignmentfile is None:
        raise ValueError('Supply alignment (BAM) files')

    if args.featureTags is None:
        raise ValueError('Supply features to extract. Supply --showtags to get a list of available tags')

    featureTags = args.featureTags.split(',')
    countTable = collections.defaultdict(
        collections.Counter)  # cell->feature->count

    def tagToHumanName(tag, TagDefinitions):
        if tag not in TagDefinitions:
            return tag
        return TagDefinitions[tag].humanName

    if args.o is not None:
        o = args.o
        tf = open(args.o, 'w')
        # Write header:
        tf.write('\t'.join([tagToHumanName(t, TagDefinitions)
                            for t in featureTags]) + '\n')
    wrote = 0
    try:
        with pysam.AlignmentFile(args.alignmentfile, ignore_truncation=True) as f:
            for i, read in enumerate(f):
                if args.dedup and read.is_duplicate:
                    continue
                if args.noqcfail and read.is_qcfail:
                    continue

                values = [
                    str(singlecellmultiomics.modularDemultiplexer.metaFromRead(read, tag))
                    for tag in featureTags
                ]
                if args.e and any((v is None for v in values)):
                    continue
                line = '%s\n' % '\t'.join(values)
                if args.o is None:
                    print(line, end="")
                else:
                    tf.write(line)
                wrote += 1
                if args.head and wrote > args.head:
                    break
    except (KeyboardInterrupt, BrokenPipeError) as e:
        pass

    if args.o is not None:
        tf.close()
