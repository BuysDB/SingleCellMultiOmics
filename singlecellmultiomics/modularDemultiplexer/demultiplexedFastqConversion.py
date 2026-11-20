#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import collections
import argparse
import singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods
from singlecellmultiomics.fastqProcessing.fastqIterator import FastqIterator


def getFeatureValue(tag, record, tags):
    """ Get the value of a requested sequence/quality feature for one record.
    SEQ / QUAL refer to the record itself, any other feature refers to a tag in the demultiplexed header.
    """
    if tag == 'SEQ':
        return record.sequence
    if tag == 'QUAL':
        return record.qual
    value = tags.get(tag)
    if value is None:
        raise ValueError(
            f'Requested feature {tag} is not present in the header of read: '
            f'{record.header}. Tags found in this record: '
            f'{",".join(sorted(tags))}')
    return value

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract features from a demultiplexed fastq file.
     Example usage: -seqfeatures BC,RX,SEQ -qualityfeatures QT,RQ,QUAL | gzip > bc_umi.fastq.gz

     """)  # -o ./bc_umi.fastq.gz
    argparser.add_argument(
        '-seqfeatures',
        type=str,
        help="Features to put into the ",
        required=True)
    argparser.add_argument(
        '-qualityfeatures',
        type=str,
        help="Features to put into ",
        required=True)
    argparser.add_argument(
        'fastqfiles',
        nargs='*',
        type=str,
        help="Input demultiplexed fastq files")

    argparser.add_argument('-format', type=str, help="Output format: fq, tsv",
                           required=True, choices=['fq', 'tsv'])
    #argparser.add_argument('-o',  type=str, help="output file path", required=True)
    args = argparser.parse_args()

    seqfeatures = args.seqfeatures.split(',')
    qualityfeatures = args.qualityfeatures.split(',')
    if len(seqfeatures) != len(qualityfeatures):
        raise ValueError('Supply a quality feature for every sequence feature')

    TagDefinitions = singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods.TagDefinitions
    for reads in FastqIterator(*args.fastqfiles):
        for readIndex, record in enumerate(reads):
            tr = singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods.TaggedRecord(
                TagDefinitions)
            tr.fromTaggedFastq(record)
            sequence = ''.join(
                (getFeatureValue(tag, record, tr.tags) for tag in seqfeatures))
            qualities = ''.join(
                (getFeatureValue(tag, record, tr.tags) for tag in qualityfeatures))
            if len(sequence) != len(qualities):
                raise ValueError(
                    'Length of base qualities does not match the length of the selected sequence features')

            if args.format == 'fq':
                print(f'{record.header}\n{sequence}\n+\n{qualities}')
            elif args.format == 'tsv':
                raise NotImplementedError(
                    'tsv output is not implemented yet, use -format fq')
