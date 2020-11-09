#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import singlecellmultiomics.modularDemultiplexer
import os
import sys
import pysam
import collections
import argparse
import pandas as pd
import numpy as np
import itertools
import singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods
import gzip  # for loading blacklist bedfiles
TagDefinitions = singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods.TagDefinitions


def coordinate_to_sliding_bin_locations(dp, bin_size, sliding_increment):
    """
    Convert a single value to a list of overlapping bins

    Parameters
    ----------
    point : int
        coordinate to look up

    bin_size : int
        bin size

    sliding_increment : int
        sliding window offset, this is the increment between bins

    Returns
    -------
    start : int
        the start coordinate of the first overlapping bin
    end :int
        the end of the last overlapping bin

    start_id : int
        the index of the first overlapping bin
    end_id : int
        the index of the last overlapping bin

    """
    start_id = int(np.ceil(((dp - bin_size) / sliding_increment)))
    start = start_id * sliding_increment
    end_id = int(np.floor(dp / sliding_increment))
    end = end_id * sliding_increment + bin_size
    return start, end, start_id, end_id


def coordinate_to_bins(point, bin_size, sliding_increment):
    """
    Convert a single value to a list of overlapping bins

    Parameters
    ----------
    point : int
        coordinate to look up

    bin_size : int
        bin size

    sliding_increment : int
        sliding window offset, this is the increment between bins

    Returns
    -------
    list: [(bin_start,bin_end), .. ]

    """
    start, end, start_id, end_id = coordinate_to_sliding_bin_locations(
        point, bin_size, sliding_increment)
    return [(i * sliding_increment, i * sliding_increment + bin_size)
            for i in range(start_id, end_id + 1)]


def read_has_alternative_hits_to_non_alts(read):
    if read.has_tag('XA'):
        for alt_align in read.get_tag('XA').split(';'):
            if len(alt_align) == 0:  # Sometimes this tag is empty for some reason
                continue

            hchrom, hpos, hcigar, hflag = alt_align.split(',')
            if not hchrom.endswith('_alt'):
                return True
    return False


def readTag(read, tag, defective='None'):
    try:
        value = singlecellmultiomics.modularDemultiplexer.metaFromRead(
            read, tag)
    except Exception as e:
        value = defective
    return value

# define if a read should be used


def read_should_be_counted(read, args, blacklist_dic = None):
    """
    Check if a read should be counted given the filter arguments

    Parameters
    ----------
    read : pysam.AlignedSegment or None
        read to check if it should be counted

    Returns
    -------
    bool
    """

    if args.r1only and read.is_read2:
        return False
    if args.r2only and read.is_read1:
        return False

    if args.filterMP:
        if not read.has_tag('mp'):
            return False

        if read.get_tag('mp')!='unique':
            return False

    if read is None or read.is_qcfail:
        return False

    # Mapping quality below threshold
    if read.mapping_quality < args.minMQ:
        return False


    if args.proper_pairs_only and not read.is_proper_pair:
        return False

    if args.no_indels and ('I' in read.cigarstring or 'D' in read.cigarstring):
        return False

    if args.max_base_edits is not None and read.has_tag('NM') and int(read.get_tag('NM'))>args.max_base_edits:
        return False

    if args.no_softclips and 'S' in read.cigarstring:
        return False

    # Read has alternative hits
    if args.filterXA:
        if read_has_alternative_hits_to_non_alts(read):
            return False

    # Read is a duplicate
    # (args.dedup and ( not read.has_tag('RC') or (read.has_tag('RC') and read.get_tag('RC')!=1))):
    if read.is_unmapped or (args.dedup and (
            read.has_tag("RR") or read.is_duplicate)):
        return False

    # Read is in blacklist
    if blacklist_dic is not None:
        if read.reference_name in blacklist_dic:
            # iterate through tuples of startend to check
            for startend in blacklist_dic[read.reference_name]:
                # start is 0-based inclusive, end is 0-based exclusive
                start_bad = read.reference_start >= startend[0] and read.reference_start < startend[1]
                end_bad = read.reference_end >= startend[0] and read.reference_end < startend[1]
                if start_bad or end_bad:
                    return False

    return True


def tagToHumanName(tag, TagDefinitions):
    if tag not in TagDefinitions:
        return tag
    return TagDefinitions[tag].humanName


def assignReads(
        read,
        countTable,
        args,
        joinFeatures,
        featureTags,
        sampleTags,
        more_args=[],
        blacklist_dic = None):

    assigned = 0
    if not read_should_be_counted(read, args, blacklist_dic):
        return assigned

    # Get the sample to which this read belongs
    sample = tuple(readTag(read, tag) for tag in sampleTags)

    # Decide how many counts this read yields
    countToAdd = 1

    if args.r1only or args.r2only:
        countToAdd = 1
    elif not args.doNotDivideFragments: # not not = True
        # IF the read is paired, and the mate mapped, we should count 0.5, and will end up
        # with 1 in total

        countToAdd = (0.5 if (read.is_paired and not read.mate_is_unmapped) else 1)

    assigned += 1

    if args.divideMultimapping:
        if read.has_tag('XA'):
            countToAdd = countToAdd / len(read.get_tag('XA').split(';'))
        elif read.has_tag('NH'):
            countToAdd = countToAdd / int(read.get_tag('NH'))
        else:
            countToAdd = countToAdd

    # Define what counts to add to what samples
    # [ (sample, feature, increment), .. ]
    count_increment = []

    feature_dict = {}
    joined_feature = []
    used_features = []
    for tag in featureTags:
        feat = str(readTag(read, tag))
        feature_dict[tag] = feat
        if args.bin is not None and args.binTag == tag:
            continue
        if args.byValue is not None and tag == args.byValue:
            continue
        joined_feature.append(feat)
        used_features.append(tag)

    if joinFeatures:
        if args.splitFeatures:
            # Obtain all key states:
            key_states = []
            for tag in featureTags:
                value = feature_dict.get(tag)
                key_states.append(value.split(args.featureDelimiter))

            for state in itertools.product(*key_states):
                joined_feature = []
                feature_dict = {
                    feature: value
                    for feature, value in zip(featureTags, state)
                }
                for feature, value in zip(featureTags, state):
                    if args.bin is not None and args.binTag == feature:
                        continue

                    if len(value) > 0:
                        joined_feature.append(value)
                    else:
                        joined_feature.append('None')

                if args.byValue is not None:
                    raise NotImplementedError(
                        'By value is not implemented for --splitFeatures')

                else:
                    count_increment.append({
                        'key': tuple(joined_feature),
                        'features': feature_dict,
                        'samples': [sample],
                        'increment': countToAdd})
                    joined_feature[0]
        else:
            if args.byValue is not None:

                try:
                    add = float(feature_dict.get(args.byValue, 0))
                except ValueError:
                    add = 0

                count_increment.append({
                    'key': tuple(joined_feature),
                    'features': feature_dict,
                    'samples': [sample],
                    'increment': add})

            else:
                count_increment.append({
                    'key': tuple(joined_feature),
                    'features': feature_dict,
                    'samples': [sample],
                    'increment': countToAdd})
    else:
        if args.bin is not None:
            raise NotImplementedError('Try using -joinedFeatureTags')

        for feature, value in feature_dict.items():
            if args.byValue is not None and feature == args.byValue:

                try:
                    add = float(feature_dict.get(args.byValue, 0))
                except ValueError:
                    add = 0

                count_increment.append({
                    'key': tuple(joined_feature),
                    'features': feature_dict,
                    'samples': [sample],
                    'increment': add})

            elif args.splitFeatures:
                for f in value.split(args.featureDelimiter):
                    count_increment.append({
                        'key': (f),
                        'features': {feature: f},
                        'samples': [sample],
                        'increment': countToAdd
                    })

            else:
                count_increment.append({
                    'key': (value, ),
                    'features': {feature: value},
                    'samples': [sample],
                    'increment': countToAdd})

    """
    Now we have a list of dicts:
    {
    'key':feature,
    'features':{feature:value},
    'samples':[sample],
    'increment':countToAdd})
    }
    """

    # increment the count table accordingly:
    if args.bin is not None:
        for dtable in count_increment:
            key = dtable['key']
            countToAdd = dtable['increment']
            samples = dtable['samples']
            value_to_be_binned = dtable['features'].get(args.binTag, None)

            if value_to_be_binned is None or value_to_be_binned == 'None':
                continue

            for start, end in coordinate_to_bins(
                    int(value_to_be_binned), args.bin, args.sliding):
                # Reject bins outside boundary
                if not args.keepOverBounds and (
                        start < 0 or end > args.ref_lengths[read.reference_name]):
                    continue
                for sample in samples:
                    countTable[sample][tuple(
                        list(key) + [start, end])] += countToAdd

    elif args.bedfile is not None:

        for dtable in count_increment:
            if args.byValue:
                key = [args.byValue]
            else:
                key = dtable['key']
            countToAdd = dtable['increment']
            samples = dtable['samples']

            # Get features from bedfile
            start, end, bname = more_args[0], more_args[1], more_args[2]
            jfeat = tuple(list(key) + [start, end, bname])
            if len(key):
                countTable[sample][jfeat] += countToAdd
            # else: this will also emit non assigned reads
            #    countTable[sample][ 'None' ] += countToAdd

    else:
        for dtable in count_increment:
            key = dtable['key']
            countToAdd = dtable['increment']
            samples = dtable['samples']

            for sample in samples:
                if len(key) == 1:
                    countTable[sample][key[0]] += countToAdd
                else:
                    countTable[sample][key] += countToAdd

    return assigned


def create_count_table(args, return_df=False):

    if len(args.alignmentfiles) == 0:
        raise ValueError("Supply at least one bam file")
    if args.bedfile is not None:
        assert(os.path.isfile(args.bedfile))

    if args.sliding is None:
        args.sliding = args.bin

    if not return_df and args.o is None and args.alignmentfiles is not None:
        args.showtags = True

    if args.showtags:
        # Find which tags are available in the file:
        head = 1000
        tagObs = collections.Counter()
        for bamFile in args.alignmentfiles:
            with pysam.AlignmentFile(bamFile) as f:
                for i, read in enumerate(f):
                    tagObs += collections.Counter([k for k,
                                                   v in read.get_tags(with_value_type=False)])
                    if i == (head - 1):
                        break
        import colorama

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

        print(f'{colorama.Style.BRIGHT}\nAVO lab tag list:{colorama.Style.RESET_ALL}')
        for tag, t in TagDefinitions.items():
            print(f'{colorama.Style.BRIGHT}{tag}{colorama.Style.RESET_ALL}\t{t.humanName}\t{colorama.Style.DIM}{"PHRED" if t.isPhred else ""}{colorama.Style.RESET_ALL}')
        exit()

    if args.o is None and not return_df:
        raise ValueError('Supply an output file')

    if args.alignmentfiles is None:
        raise ValueError('Supply alignment (BAM) files')

    joinFeatures = False
    featureTags = []
    if args.featureTags is not None:
        featureTags = args.featureTags.split(',')

    if args.joinedFeatureTags is not None:
        featureTags = args.joinedFeatureTags.split(',')
        joinFeatures = True

        # When -byValue is used, and joined feature tags are supplied, automatically append the args.byValue tag, otherwise it will get lost
        if args.byValue is not None and len(featureTags)>0 and args.byValue not in featureTags:
            featureTags.append(args.byValue)


    if args.bin is not None and args.binTag not in featureTags:
        print("The bin tag was not supplied as feature, automatically appending the bin feature.")
        featureTags.append(args.binTag)

    if len(featureTags) == 0:
        raise ValueError(
            "No features supplied! Please supply -featureTags -joinedFeatureTags and or -binTag")

    sampleTags = args.sampleTags.split(',')
    countTable = collections.defaultdict(
        collections.Counter)  # cell->feature->count

    if args.blacklist is not None:
        # create blacklist dictionary {chromosome : [ (start1, end1), ..., (startN, endN) ]}
        # used to check each read and exclude if it is within any of these start end sites
        #
        blacklist_dic = {}
        print("Creating blacklist dictionary:")
        with open(args.blacklist, mode='r') as blfile:
            for row in blfile:
                parts = row.strip().split()
                chromo, start, end = parts[0], int(parts[1]), int(parts[2])
                if chromo not in blacklist_dic:
                    # init chromo
                    blacklist_dic[chromo] = []  # list of tuples
                blacklist_dic[chromo].append( (start, end) )
        print(blacklist_dic)
    else:
        blacklist_dic = None

    assigned = 0
    for bamFile in args.alignmentfiles:

        with pysam.AlignmentFile(bamFile) as f:
            i = 0  # make sure i is defined
            if args.bin:
                # Obtain the reference sequence lengths
                ref_lengths = {
                    r: f.get_reference_length(r) for r in f.references}
                args.ref_lengths = ref_lengths
            if args.bedfile is None:
                # for adding counts associated with a tag OR with binning
                if args.contig is not None:
                    pysam_iterator = f.fetch(args.contig)
                else:
                    pysam_iterator = f

                for i, read in enumerate(pysam_iterator):
                    if i % 1_000_000 == 0:
                        print(
                            f"{bamFile} Processed {i} reads, assigned {assigned}, completion:{100*(i/(0.001+f.mapped+f.unmapped+f.nocoordinate))}%")

                    if args.head is not None and i > args.head:
                        break

                    assigned += assignReads(read,
                                            countTable,
                                            args,
                                            joinFeatures,
                                            featureTags,
                                            sampleTags,
                                            blacklist_dic = blacklist_dic)
            else:  # args.bedfile is not None
                # for adding counts associated with a bedfile
                with open(args.bedfile, "r") as bfile:
                    #breader = csv.reader(bfile, delimiter = "\t")
                    for row in bfile:

                        parts = row.strip().split()
                        chromo, start, end, bname = parts[0], int(float(
                            parts[1])), int(float(parts[2])), parts[3]
                        if args.contig is not None and chromo != args.contig:
                            continue
                        for i, read in enumerate(f.fetch(chromo, start, end)):
                            if i % 1_000_000 == 0:
                                print(
                                    f"{bamFile} Processed {i} reads, assigned {assigned}, completion:{100*(i/(0.001+f.mapped+f.unmapped+f.nocoordinate))}%")
                            assigned += assignReads(read,
                                                    countTable,
                                                    args,
                                                    joinFeatures,
                                                    featureTags,
                                                    sampleTags,
                                                    more_args=[start,
                                                               end,
                                                               bname],
                                                    blacklist_dic = blacklist_dic)

                            if args.head is not None and i > args.head:
                                break

            print(
                f"Finished: {bamFile} Processed {i} reads, assigned {assigned}")
    print(f"Finished counting, now exporting to {args.o}")
    df = pd.DataFrame.from_dict(countTable)

    # Set names of indices
    if not args.noNames:
        df.columns.set_names([tagToHumanName(t, TagDefinitions)
                              for t in sampleTags], inplace=True)

        try:
            if args.bin is not None:
                index_names = [tagToHumanName(
                    t, TagDefinitions) for t in featureTags if t != args.binTag] + ['start', 'end']
                df.index.set_names(index_names, inplace=True)
            elif args.bedfile is not None:
                index_names = [tagToHumanName(
                    t, TagDefinitions) for t in featureTags if t != args.binTag] + ['start', 'end', 'bname']
                df.index.set_names(index_names, inplace=True)
            elif joinFeatures:
                index_names = [
                    tagToHumanName(
                        t, TagDefinitions) for t in featureTags]
                df.index.set_names(index_names, inplace=True)
            else:
                index_names = ','.join(
                    [tagToHumanName(t, TagDefinitions) for t in featureTags])
                df.index.set_names(index_names, inplace=True)
        except Exception as e:
            pass
        print(index_names)

    if return_df:
        return df

    if args.o.endswith('.pickle') or args.o.endswith('.pickle.gz'):
        df.to_pickle(args.o)
    else:
        df.to_csv(args.o)
    return args.o
    print("Finished export.")


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Convert BAM file(s) which has been annotated using featureCounts or otherwise to a count matrix')
    argparser.add_argument(
        '-o',
        type=str,
        help="output csv path, or pandas dataframe if path ends with pickle.gz",
        required=False)
    argparser.add_argument(
        '-featureTags',
        type=str,
        default=None,
        help='Tag(s) used for defining the COLUMNS of the matrix. Single dimension.')
    argparser.add_argument(
        '-joinedFeatureTags',
        type=str,
        default=None,
        help='These define the COLUMNS of your matrix. For example if you want allele (DA) and restriction site (DS) use DA,DS. If you want a column containing the chromosome mapped to use "chrom" as feature. Use this argument if you want to use a multidimensional index')
    argparser.add_argument(
        '-sampleTags',
        type=str,
        default='SM',
        help='Comma separated tags defining the names of the ROWS in the output matrix')
    argparser.add_argument('alignmentfiles', type=str, nargs='*')
    argparser.add_argument(
        '-head',
        type=int,
        help='Run the algorithm only on the first N reads to check if the result looks like what you expect.')


    argparser.add_argument(
        '--r1only',
        action='store_true',
        help='Only count R1')

    argparser.add_argument(
        '--r2only',
        action='store_true',
        help='Only count R2')

    argparser.add_argument(
        '--splitFeatures',
        action='store_true',
        help='Split features by , . For example if a read has a feature Foo,Bar increase counts for both Foo and Bar')
    argparser.add_argument('-featureDelimiter', type=str, default=',')
    argparser.add_argument(
        '-contig',
        type=str,
        help='Run only on this chromosome')

    multimapping_args = argparser.add_argument_group('Multimapping', '')
    multimapping_args.add_argument(
        '--divideMultimapping',
        action='store_true',
        help='Divide multimapping reads over all targets. Requires the XA or NH tag to be set.')
    multimapping_args.add_argument(
        '--doNotDivideFragments',
        action='store_true',
        help='When used every read is counted once, a fragment will count as two reads. 0.5 otherwise')


    binning_args = argparser.add_argument_group('Binning', '')
    #binning_args.add_argument('-offset', type=int, default=0, help="Add offset to bin. If bin=1000, offset=200, f=1999 -> 1200. f=4199 -> 3200")
    binning_args.add_argument(
        '-sliding',
        type=int,
        help="make bins overlapping, the stepsize is equal to the supplied value here. If nothing is supplied this value equals the bin size")
    binning_args.add_argument(
        '--keepOverBounds',
        action='store_true',
        help="Keep bins which go over chromsome bounds (start<0) end > chromsome length")

    binning_args.add_argument(
        '-bin',
        type=int,
        help="Devide and floor to bin features. If bin=1000, f=1999 -> 1000.")
    # binning_args.add_argument('--showBinEnd', action='store_true', help="If
    # True, then show DS column as 120000-220000, otherwise 120000 only. This
    # specifies the bin range in which the read was counted" ) this is now
    # always on!
    binning_args.add_argument('-binTag', default='DS')
    binning_args.add_argument(
        '-byValue',
        type=str,
        help='Extract the value from the supplied tag and use this as count to add')
    binning_args.add_argument('-blacklist', metavar='BED FILE BLACKLIST TO FILTER OUT', help = 'Bedfile of blacklist regions to exclude', default=None)
    bed_args = argparser.add_argument_group('Bedfiles', '')
    bed_args.add_argument(
        '-bedfile',
        type=str,
        help="Bed file containing 3 columns, chromo, start, end to be read for fetching counts")



    filters = argparser.add_argument_group('Filters', '')

    filters.add_argument(
        '--proper_pairs_only',
        action='store_true',
        help='Only count reads mapped in a proper pair (within the expected insert size range)')

    filters.add_argument(
        '--no_softclips',
        action='store_true',
        help='Only count reads without softclips')

    filters.add_argument(
        '-max_base_edits',
        type=int,
        help='Count reads with at most this value of bases being different than the reference')


    filters.add_argument(
        '--no_indels',
        action='store_true',
        help='Only count reads without indels')

    filters.add_argument(
        '--dedup',
        action='store_true',
        help='Count only the first occurence of a molecule. Requires RC tag to be set. Reads without RC tag will be ignored!')

    filters.add_argument(
        '-minMQ',
        type=int,
        default=0,
        help="minimum mapping quality")

    filters.add_argument(
        '--filterMP',
        action= 'store_true',
        help="Filter reads which are not uniquely mappable, this is based on presence on the `mp` tag, which can be generated by bamtagmultiome.py ")

    filters.add_argument(
        '--filterXA',
        action='store_true',
        help="Do not count reads where the XA (alternative hits) tag has been set for a non-alternative locus.")

    argparser.add_argument(
        '--noNames',
        action='store_true',
        help='Do not set names of the index and columns, this simplifies the resulting CSV/pickle file')
    argparser.add_argument(
        '--showtags',
        action='store_true',
        help='Show a list of commonly used tags')

    args = argparser.parse_args()

    create_count_table(args)
