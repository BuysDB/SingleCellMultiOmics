#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
import collections
import pysamiterators.iterators as pysamIterators

complement = str.maketrans('ATCGN', 'TAGCN')

"""
This module contains functions which use data encoded in bam or fastq tags
"""


# For a mapped read pair it is very important to figure out which bases are actual genomic signal
# Al sequence behind and on top of the random primer cannot be trusted and should be masked out
# secondly all signal before the starting location of R1 cannot be trusted
# Due to dovetailing this cannot be detected before mapping
# This function returns a lower and higher bound of the locations within the fragment that can be trusted
# ASCII art: (H is primer sequence)
#           R1 H------------------------>
#   <------------------HH R2

# Output: (E is emitted)
#           R1 HEEEEEEE----------------->
#   <------------------HH R2

def getPairGenomicLocations(R1, R2, R1PrimerLength=4, R2PrimerLength=6):

    if R1 is None or R2 is None or R1.is_unmapped:
        # This is an annoying situation, we cannot determine what bases can be
        # trusted
        raise ValueError('Genomic locations cannot be determined')
    if R2.is_unmapped:
        # This is an annoying situation, we cannot determine what bases can be
        # trusted
        raise ValueError('Genomic locations cannot be determined')

    if R1.is_reverse == R2.is_reverse:
        # The fragment is not correctly mapped
        raise ValueError('Fragment incorrectly mapped')

    if not R1.is_reverse:  # R1 is forward, R2 is reverse
        #           R1 H------------------------>
        #   <------------------HH R2
        start = R1.reference_start + R1PrimerLength
        end = R2.reference_end - R2PrimerLength
    else:
        #           R2 HH------------------------>
        #   <------------------HH R1
        start = R2.reference_start + R2PrimerLength
        end = R1.reference_end - R1PrimerLength

    if start >= end:
        raise ValueError('Fragment has no size')

    return start, end


# Get an identifier which identifies the reverse transcription reaction
# (hexamerStart, hexamerBases)
def getRandomPrimerHash(R2, onStart=True, primerLength=6):
    if R2.query_sequence is None:
        return None, None
    # The read was not mapped
    if R2.is_unmapped:
        # Guess the orientation does not matter
        return None, R2.query_sequence[:primerLength]

    if R2.is_reverse:
        global complement
        return(R2.reference_end, R2.query_sequence[-primerLength:][::-1].translate(complement))
    else:
        return(R2.reference_start, R2.query_sequence[:primerLength])

    raise ValueError()


def getListAllPositions(l):
    positions = set()

    for read in l:
        if isinstance(read, list):
            for x in read:
                for start, end in x.get_blocks():
                    for x in range(start, end):
                        positions.add(x)
        else:
            for start, end in read.get_blocks():
                for x in range(start, end):
                    positions.add(x)
    return positions


def getCoverageAllPositions(l):
    positions = collections.Counter()
    for read in l:
        for start, end in read.get_blocks():
            for x in range(start, end):
                positions[x] += 1
    return positions


def getMateDictSpanningCoordinates(md):
    surfaceStart = None
    surfaceEnd = None
    contig = None
    for fragmentId, reads in md.items():
        for read in reads.values():
            if contig is None and read.reference_name is not None:
                contig = read.reference_name
            if surfaceStart is None or read.reference_start < surfaceStart:
                surfaceStart = read.reference_start
            if surfaceEnd is None or read.reference_end > surfaceEnd:
                surfaceEnd = read.reference_end
    return contig, surfaceStart, surfaceEnd


def getUniqueRandomPrimers(readIter, primerLength=6):
    randomPrimers = set()
    for read in readIter:
        if read.is_read2:
            randomPrimers.add(
                getRandomPrimerHash(
                    read, primerLength=primerLength))
    return randomPrimers

# Obtain from how many sources the reads are derived
# The sequencer, lane and tile are taken into account for this metric


def getSources(readIter):
    instruments = set()
    lanes = set()  # instrument, flowcell, lane
    tiles = set()  # instrument, lane, tile
    for read in readIter:
        instrument = read.get_tag('Is')
        flowcell = read.get_tag('Fc')
        lane = read.get_tag('La')
        tile = read.get_tag('Ti')
        instruments.add(instrument)
        lanes.add((instrument, flowcell, lane))
        tiles.add((instrument, flowcell, lane, tile))
    return instruments, lanes, tiles


def getCycleOffset(read):

    start = 0
    if read.is_read1:
        start = read.get_tag('eB') if read.has_tag('eB') else 0
    elif read.is_read2:
        start = read.get_tag('EB') if read.has_tag('EB') else 0
    else:
        raise ValueError('Designed for single or mate pair only')
    return (start)


def getReadTotalCycles(read, cycleOffset=None):
    # The obvious part:
    totalCycles = read.infer_read_length()
    # Add trimmed cycles:
    if cycleOffset is None:
        cycleOffset = getCycleOffset(read)
    totalCycles += cycleOffset  # @warn: end is not defined!
    return totalCycles

# This iterator is similar and a wrapper of the Pysam get_aligned_pairs function
# The difference is that the cycle of the sequencer is emitted (distToFirstCycle) (int or float)
# yields cycle, queryIndex, referencePos, (refbase)
# The cycle of the sequencer is obtained by the index, the read orientation and the eB/EB tags
# The second added feature is that a reference handle can be added which
# will yield reference bases from the supplied fasta file. This feature is
# neccesary when mapping to masked genomes


class ReadCycleIterator():

    def __init__(self,
                 read,
                 matches_only=False,  # transfered to pysam api
                 with_seq=False,  # emit reference bases
                 emitFloats=False,  # Emit as percentage of total cycles instead of absolute cycles
                 reference=None,
                 # obtain reference base from a reference (Should be type pysam
                 # FastaFile)

                 ):
        self.read = read
        self.with_seq = with_seq
        self.matches_only = matches_only
        self.emitFloats = emitFloats
        self.start = getCycleOffset(read)
        self.len = getReadTotalCycles(read, cycleOffset=self.start)
        self.reference = reference

    def __repr__(self):
        return(f'{self.read}, {self.len} cycles starting at {self.start}')

    def __iter__(self):
        if self.reference is None:
            self.iterator = iter(
                self.read.get_aligned_pairs(
                    matches_only=self.matches_only,
                    with_seq=self.with_seq))
        else:
            self.iterator = iter(
                pysamIterators.ReferenceBackedGetAlignedPairs(
                    self.read,
                    self.reference,
                    matches_only=self.matches_only,
                    with_seq=True))
        return self

    def __next__(self):

        if self.with_seq:
            readIndex, referencePos, referenceBase = next(self.iterator)
        else:
            readIndex, referencePos = next(self.iterator)
        # Obtain cycle:
        if not self.read.is_reverse:
            cycle = readIndex + self.start
        else:
            cycle = self.len - readIndex - self.start - \
                1  # minus one as the index starts at 0
        if self.emitFloats:
            if self.len == 0:
                cycle = 0
            else:
                cycle /= self.len
        if self.with_seq:
            return cycle, readIndex, referencePos, referenceBase
        else:
            return cycle, readIndex, referencePos
