#!/usr/bin/env python3
# -*- coding: utf-8 -*-
c = 1_000 # !!! PLEASE USE PYTHON 3.6 OR HIGHER !!!
import os,pysam,re
import subprocess,multiprocessing
import glob
import sys,time
import collections
import uuid
from singlecellmultiomics.alleleTools import alleleTools
import pickle
import gzip
complement = str.maketrans('ATGC', 'TACG')
import pysamiterators.iterators as pysamIterators
import argparse
from singlecellmultiomics.tagtools import tagtools
import singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods
import singlecellmultiomics.features
import math
import numpy as np

def phred_to_prob(phred):
    try:
        return math.pow(10,-(ord(phred)-33)/10 )
    except ValueError:
        return 1

def hamming_distance(a,b):
    return sum((i!=j and i!='N' and j!='N' for i,j in zip(a,b)))

if __name__ == "__main__" :
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""Add allelic NLA and/or MSPJI digest information to a demultiplexed bam file
    Adds the following tags:
    RC: number of oversequenced molecules
    DT: Source type (NLA/MSPJI/RNA.. )
    DA: Digest allele
    RP: Recognized sequence
    RS: Recognized strand '-','?','+'

    Which tags are filled in depends on which files you supply, for example if you do not supply an allelic VCF file, then no allele tag will be added
    """)

    argparser.add_argument('bamfiles', metavar='bamfile', type=str, nargs='+')
    argparser.add_argument('--local', action='store_true', help='Do not qsub sorting and indexing [deprecated]')
    #argparser.add_argument('--noSort', action='store_true', help='Do not sort and index')
    #argparser.add_argument('--cluster', action='store_true', help='Do not qsub sorting and indexing')


    tagAlgs = argparser.add_argument_group('Tagging algorithms', '')
    tagAlgs.add_argument('--ftag', action='store_true', help='Add query name tags added by demultiplexer')
    tagAlgs.add_argument('--atag', action='store_true', help='Add only allele tags')
    tagAlgs.add_argument('--mspji', action='store_true', help='Look for mspji digest')
    tagAlgs.add_argument('--nla', action='store_true', help='Look for nlaIII digest')
    tagAlgs.add_argument('--chic', action='store_true', help='Look for mnase (scCHIC) digest')
    tagAlgs.add_argument('--scar', action='store_true', help='Create R1 start, cigar sequence based DS tags')
    tagAlgs.add_argument('-tag', type=str, default=None, help='Determine oversequencing based on a tag (for example XT to count RNA oversequencing for featureCounts counted transcripts. chrom for chromosome/contig count)')


    ### RNA options
    tagAlgs.add_argument('--rna', action='store_true', help='Assign RNA molecules, requires a intronic and exonic GTF file')
    argparser.add_argument('-introns', type=str, help='Intronic GTF file')
    argparser.add_argument('-exons', type=str, help='Exonic GTF file.')
    """
    How much bases of the fragment fall inside the intron(s)
    How much bases of the fragment fall inside the exon(s)
    is the fragment spliced? How many bases are spliced out?
    """

    argparser.add_argument('-moleculeRadius', type=int, help='Search radius for each molecule, if a cut site is found within this range with the same umi, it will be counted as the same molecule', default=0)

    argparser.add_argument('-alleles', type=str, help='VCF file. Add allelic info flag based on vcf file')
    argparser.add_argument('--loadAllelesToMem',  action='store_true',help='Load allele data completely into memory')

    argparser.add_argument('-ref', type=str, help='reference FASTA file.')
    argparser.add_argument('-o', type=str, default='./tagged', help='output folder')
    argparser.add_argument('--dedup', action='store_true', help='Create deduplicated bam files')



    argparser.add_argument('-chr', help='run only on this chromosome')

    argparser.add_argument('-head', type=int, default=None, help='Only process first n reads of every bam file')

    devArgs = argparser.add_argument_group('Development options', '')
    devArgs.add_argument('--fatal', action='store_true', help='Crash upon any encountered error')
    devArgs.add_argument('--verbose', action='store_true', help='Be verbose')
    devArgs.add_argument('--knh', action='store_true', help='Keep non-headered bam file')

    args = argparser.parse_args()

    if not args.mspji and not args.nla and not args.chic and not args.ftag  and not args.rna and args.tag is None and args.atag is None:
        raise ValueError('Please supply any or a combination of --ftag --nla --chic --mspji')

    if args.rna and (args.introns is None or args.exons is None):
        raise ValueError('Please supply exonic and intronic GTF files -introns -exons')

    """
    Corrected Tags
    BC: This should be used for the sample-barcode so that there is no conflict in case that both sample- and molecular-barcodes are used
    QT: This should be the quality for the sample-barcode BC

    New Tags
    RX: The raw sequence bases of the molecular-barcode(s) of the read(s) can be placed in this tag
    QX: The original qualites (as QUAL) of the bases in RX

    BX: Sequence bases of the unique molecular identifier, may be corrected or raw
    BZ: Quality score the unique molecular identifier. May be corrected or raw. See BX tag for qualities

    MI: Molecular-barcode sequence or identifier, may include non-sequence bases (Here, it is BC+RX)
    QM: Molecular-barcode qualities, QT+QX

    DS: coordinate of site
    RC: number of oversequenced molecules
    DT: Source type (NLA/MSPJI/RNA..? )
    DA: Digest allele
    RP: Recognized sequence
    RS: Recognized strand '-','?','+'
    LI: ligation sequence

    """


class RNA_Flagger():

    def __init__(self, reference=None, alleleResolver=None, moleculeRadius=0, verbose=False, exon_gtf=None, intron_gtf=None, **kwargs):


        self.annotations= {}
        self.annotations['EX'] = singlecellmultiomics.features.FeatureContainer()
        self.annotations['EX'].loadGTF( exon_gtf, select_feature_type=['exon'] )
        self.annotations['IN'] = singlecellmultiomics.features.FeatureContainer()
        self.annotations['IN'].loadGTF( intron_gtf, select_feature_type=['intron'] )

        self.exon_hit_tag = 'eH'
        self.intron_hit_tag = 'iH'
        self.assinged_exons_tag = 'AE'
        self.assinged_introns_tag = 'AI'
        self.is_spliced_tag = 'SP' #unused

        self.overlap_tag = 'XM'

    def digest(self, reads):

        feature_overlap = collections.Counter() #feature->overlap
        # Go over reads and calculate overlap with features

        exon_count = 0
        intron_count = 0

        for read in reads:
            if read is None:
                continue

            states = ['.']*read.query_length
            for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False):



                overlaps_with_intron = False
                overlaps_with_exon = False
                exon_hits = set()

                for hit in self.annotations['EX'].findFeaturesAt(chromosome=read.reference_name,lookupCoordinate=ref_pos,strand=None):
                    hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                    overlaps_with_exon = True
                    exon_hits.add(hit_id)

                intron_hits = set()
                for hit in self.annotations['IN'].findFeaturesAt(chromosome=read.reference_name,lookupCoordinate=ref_pos,strand=None):
                    hit_start, hit_end, hit_id, hit_strand, hit_ids = hit
                    overlaps_with_intron = True
                    intron_hits.add(hit_id)

                if overlaps_with_exon and not overlaps_with_intron:
                    exon_count+=1
                    if self.overlap_tag is not None:
                        states[q_pos] = 'Z'
                elif not overlaps_with_exon and overlaps_with_intron:
                    intron_count+=1
                    if self.overlap_tag is not None:
                        states[q_pos] = 'X'

            if self.overlap_tag is not None:
                read.set_tag(self.overlap_tag, ''.join(states))

        for read in reads:
            if read is not None:
                read.set_tag(self.exon_hit_tag, exon_count)
                read.set_tag(self.intron_hit_tag, intron_count)

                read.set_tag( self.assinged_exons_tag, ','.join( (str(e) for e in exon_hits) ))
                read.set_tag( self.assinged_introns_tag, ','.join( (str(e) for e in intron_hits) ))





def molecule_to_random_primer_dict(molecule, primer_length=6, primer_read=2): #1: read1 2: read2
    rp = collections.defaultdict(list)
    for fragment in molecule:
        if fragment[primer_read-1] is not None:
            hstart, hseq = tagtools.getRandomPrimerHash(fragment[primer_read-1], onStart=True, primerLength=6)
            rp[hstart, hseq].append(fragment)
    return rp
    #for k,p in rp.items():
    #    yield p


class MoleculeIterator():

    def __init__(self,
        alignmentfile,
        look_around_radius=100_000,
        umi_hamming_distance=0, # 0: only exact match, 1: single base distance
        sample_select=None, # iterable of samples to only select molecules from
        # when a string is supplied only a single sample is selected

        **pysam_kwargs
        ):
        self.alignmentfile = alignmentfile
        self.look_around_radius = look_around_radius

        self.current_position = None #
        self.current_chromosome = None #
        self.molecules_yielded = 0
        self.umi_hamming_distance = umi_hamming_distance
        self.pysam_kwargs = pysam_kwargs
        self.sample_select = sample_select
        self.sample_tag = 'SM'
        self._clear()

        self.filter_function = None # function which results given two reads if it is usable


    def _clear(self):
        self.molecule_cache = collections.defaultdict(
            lambda: collections.defaultdict(list)) # position -> cell,umi,strand,allele.. -> reads

    def sample_assignment_function(self, fragment):
        for read in fragment:
            if read is not None:
                if read.has_tag(self.sample_tag):
                    return read.get_tag(self.sample_tag)
        return None

    # Returns postion unique identifier for a fragment
    def localisation_function(self, fragment):
        if not fragment[0].has_tag('DS'):
            return None
        return fragment[0].get_tag('DS')


    def __repr__(self):
        return f"""Molecule iterator
        current position: {self.current_chromosome}:{self.current_position}
        yielded {self.molecules_yielded} so far
        currently {len(self.molecule_cache)} positions cached
        keeping {self.get_cached_fragment_count()} fragments cached
        """

    def get_cached_fragment_count(self):
        total_fragments = 0
        for position, data_per_molecule in self.molecule_cache.items():
            for molecule_id, molecule_reads in data_per_molecule.items():
                total_fragments+=len(molecule_reads)
        return total_fragments

    def assignment_function(self, fragment):
        return  fragment[0].get_tag('SM'),fragment[0].get_tag('RX'),fragment[0].get_tag('RS')

    # Returns True if two fragment ids  are identical or close enough
    def eq_function(self, assignment_a, assignment_b):
        sample_A, umi_A, strand_A = assignment_a
        sample_B, umi_B, strand_B = assignment_b

        if sample_A!=sample_B or strand_A!=strand_B:
            return False

        if self.umi_hamming_distance==0:
            return umi_A==umi_B
        else:
            return hamming_distance(umi_A,umi_B)<=self.umi_hamming_distance

    def _yield_all_in_current_cache(self):
        for position, data_per_molecule in self.molecule_cache.items():
            for molecule_id, molecule_reads in data_per_molecule.items():
                self.molecules_yielded +=1
                yield molecule_reads
        self._clear()

    def _purge(self):
        drop= []
        for pos in self.molecule_cache:
            if pos<(self.current_position-self.look_around_radius):
                # purge this coordinate:
                for molecule_id, molecule_reads in self.molecule_cache[pos].items():
                    self.molecules_yielded +=1
                    yield molecule_reads
                drop.append(pos)

        for pos in drop:
            del self.molecule_cache[pos]

    def __iter__(self):
        for fragment in pysamIterators.MatePairIterator( self.alignmentfile, **self.pysam_kwargs ):

            # now yield fragments which are finished :
            if fragment[0].reference_name is None:
                continue

            # Check if we want to drop the fragment because its corresponding sample is not selected
            if self.sample_select is not None:
                sample = self.sample_assignment_function(fragment)
                if type(self.sample_select)==str: # single sample
                    if self.sample_select!=sample:
                        continue
                else: # list or set of samples
                    if not sample in self.sample_select:
                        continue

            # We went to another chromsome, purge all in cache:
            if fragment[0].reference_name!=self.current_chromosome and self.current_chromosome is not None:
                for molecule in self._yield_all_in_current_cache():
                    yield molecule

            # Check if we need to purge and yield results:
            for molecule in self._purge():
                yield molecule

            position = self.localisation_function(fragment)
            if position is None:
                continue
            molecule_id = self.assignment_function(fragment)

            # Check if there is a molecule to assign to already present:
            assigned = False
            for existing_molecule in self.molecule_cache[position]:
                if self.eq_function(existing_molecule, molecule_id):
                    assigned  = True
                    self.molecule_cache[position][existing_molecule].append(fragment)
            if not assigned:
                self.molecule_cache[position][molecule_id].append(fragment)

            self.current_chromosome = fragment[0].reference_name
            self.current_position = fragment[0].reference_end

        # yield everything which was not yielded yet
        for molecule in self._yield_all_in_current_cache():
            yield molecule

class RangeCache():

    def __init__(self, maxRange=1200):
        self.d = dict() #Chrom -> pos -> value
        self.current = 0
        self.maxRange = maxRange
        self.freedMemory = 0

    # Obtain record within a radius (center, Radius) , returns first (closest) instance it finds
    def getWithinRange(self, chrom, center, radius):
        if not chrom in self.d:
            return None
        for i in range(0, radius+1):
            r = self.get(chrom, i+center)
            if r is not None:
                return r
            r = self.get(chrom, center-i)
            if r is not None:
                return r

    def get(self, chrom, pos ):
        if not chrom in self.d:
            return None

        if not pos in self.d[chrom]:
            return None
        return self.d[chrom][pos]

    def set(self, chrom, pos, data):
        self.purge(chrom, pos)
        if not chrom in self.d:
            self.d[chrom] = {}
        self.d[chrom][pos] = data

    def purge( self, chrom, pos ):
        if chrom in self.d:
            drop = []
            for position in self.d[chrom]:
                if abs(position-pos)>self.maxRange:
                    drop.append(position)
            for d in drop:
                if d in self.d[chrom]:
                    del self.d[chrom][d]
                    self.freedMemory+=1
        # Remove all data on other chromosomes:
        drop = [x for x in self.d if x!=chrom]
        for d in drop:
            del self.d[d]
            self.freedMemory+=1

    def __len__(self):
        return sum( [len(self.d[chrom]) for chrom in self.d])

# A digest function should fullfill the following conditions:
# Accept : a read pair of R1, and R2 optional
# a reference genome handle (optional)
# an allele tool handle (optional)

class DigestFlagger():

    def __init__(self, reference=None, alleleResolver=None, moleculeRadius=0, verbose=False, **kwargs):
        self.reference = reference
        self.alleleResolver = alleleResolver
        self.verbose = verbose
        self.siteCoordinateTag = 'DS'
        self.oversequencedMoleculeTag = 'RC'
        self.recognitionSequenceTag = 'RZ'
        self.sourceTypeTag = 'DT' # DataType
        self.alleleTag = 'DA'
        self.alleleSetTag = 'AA'
        self.strandTag = 'RS'
        self.umiTag = 'RX'
        self.sampleTag = 'SM'
        self.ligationTag = 'LI'
        self.rejectionTag = 'RR'
        self.randomPrimerPosTag = 'rP'
        self.randomPrimerPosSeq = 'rS'

        self.fragmentSizeTag = 'fS'
        self.fragmentEndTag = 'fe'
        self.fragmentStartTag = 'fs'

        self.cacheWrites = 0

        self.cachePurgeEvery = 10_000 # Clean the cache every N writes

        self.moleculeRadius = moleculeRadius
        # This hold which sites we saw to count oversequencing
        self.observedSites = collections.defaultdict( RangeCache ) # sample -> (chrom,allele,strand,..) -> pos -> seen

    def __repr__(self):
        return f'Tagger with {self.cacheWrites} cache writes, total consumption: {self.getTotalConsumption()} '

    def getTotalConsumption(self):
        return sum( [len(self.observedSites[sample]) for sample in self.observedSites])


    def setRejectionReason(self, read, reason):
        self.appendTag(read, self.rejectionTag, reason)

    def setLigationSite(self, read, site):
        read.set_tag( self.ligationTag, site  )

    def setSiteCoordinate( self, read, coordinate ):
        self.appendTag(read, self.siteCoordinateTag, coordinate )

    def setSiteOversequencing( self, read, moleculeIndex=1 ): # 1 if first seen 2, second, -1 if None
        read.set_tag( self.oversequencedMoleculeTag, -1 if moleculeIndex is None else moleculeIndex  )
        # Decribe as string and set tag:

        if moleculeIndex>1:
            read.is_duplicate = True

    def setFragmentSize(self, read, size):
        read.set_tag( self.fragmentSizeTag,size)

    def setFragmentTrust(self, read, start, end):
        read.set_tag( self.fragmentStartTag,start)
        read.set_tag( self.fragmentEndTag,end)


    def setAllele( self, read, allele): # 1 if first seen 2, second, -1 if None
        read.set_tag( self.alleleTag, allele)

    def setRecognizedSequence( self, read, sequence ):
        self.appendTag(read, self.recognitionSequenceTag, sequence )

    def setSource( self, read, source ):
        self.appendTag(read, self.sourceTypeTag, source )

    def setStrand( self, read, strand ):
        read.set_tag( self.strandTag, strand )


    def setRandomPrimer(self, R1,R2, hstart, hseq ):
        if R1 is not None:
            R1.set_tag(self.randomPrimerPosSeq, hseq)
            R1.set_tag(self.randomPrimerPosTag, hstart)
        if R2 is not None:
            R2.set_tag(self.randomPrimerPosSeq, hseq)
            R2.set_tag(self.randomPrimerPosTag, hstart)

    def appendTag(self, read, tag, value):
        if read is None:
            return

        if not read.has_tag(tag):
            read.set_tag(tag, value)
        else:
            read.set_tag(tag, f'{read.get_tag(tag)},{value}' )

    def addAlleleInfo(self, reads):
        allele = None if self.alleleResolver is None else self.alleleResolver.getAllele(reads)
        if self.verbose:
            print(allele,reads)
        if allele is not None and len(allele)>0:
            allele = ','.join(sorted(list(allele)))

            for read in reads:
                self.setAllele(read,allele)
        return allele

    # siteInfo describes the site as tuple: (allele, recognizedStrand, umiSequence, ..)
    # the site info can have arbitrary length, with as only requirement that it should be hashable
    def increaseAndRecordOversequencing( self, sample, chrom, pos, siteInfo=() ):
        cutSite = (chrom, pos)

        self.cacheWrites += 1
        if self.cacheWrites>self.cachePurgeEvery:
            for sample in self.observedSites:
                self.observedSites[sample].purge(*cutSite)
            self.cacheWrites=0

        if self.moleculeRadius!=0:
            current = self.observedSites[sample].getWithinRange(*cutSite, radius=self.moleculeRadius)
        else:
            current = self.observedSites[sample].get( *cutSite )

        # There are no matching molecules nearby:
        if current is None:
            current = collections.Counter()
            self.observedSites[sample].set( *cutSite, current )

        UMI = siteInfo
        current[UMI]+=1
        #print(current)
        return  current[UMI]

    ### ADD THIS FUNCTION:
    #def digest( self, reads )
    #


class QueryNameFlagger(DigestFlagger):
    def __init__(self, **kwargs):
        self.assignedReadGroups = set()
        DigestFlagger.__init__(self, **kwargs )

    def digest(self, reads):
        for read in reads:
            if read is None:
                continue
            if read.query_name.startswith('UMI'): # Old format
                import tagBamFile # this is not included anymore
                tagBamFile.recodeRead(read)
            else:
                tr = singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods.TaggedRecord(
                    singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods.TagDefinitions
                )

                tr.fromTaggedBamRecord(read)

                newHeader = tr.asIlluminaHeader()
                read.query_name = newHeader

                tr.tagPysamRead(read)
                rg = f"{read.get_tag('Fc') if read.has_tag('Fc') else 'NONE'}.{read.get_tag('La') if read.has_tag('La') else 'NONE'}.{read.get_tag('SM') if read.has_tag('SM') else 'NONE'}"
                self.assignedReadGroups.add(rg)
                # Add read group:
                read.set_tag('RG', rg)


class NlaIIIFlagger( DigestFlagger ):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )

    def addSite(self, reads, strand, restrictionChrom, restrictionPos ):

        if not reads[0].has_tag(self.sampleTag) or not reads[0].has_tag(self.umiTag):
            return

        sample = reads[0].get_tag(self.sampleTag)
        umi = reads[0].get_tag(self.umiTag)
        allele = None if not reads[0].has_tag(self.alleleTag) else reads[0].get_tag(self.alleleTag)
        siteInfo = tuple( [ x  for x in [strand, allele, umi] if x is not None])
        moleculeId = self.increaseAndRecordOversequencing(  sample,  restrictionChrom, restrictionPos, siteInfo=siteInfo)

        for read in reads:
            if read is None:
                continue
            self.setSiteOversequencing( read, moleculeId )
            self.setSiteCoordinate(read, restrictionPos)
            self.setSource(read, 'NLA'), {}
            if allele is not None:
                self.setAllele(read, allele)

            self.setStrand(read, '+' if strand==1 else ('-' if strand==0 else '?'))


    def digest(self, reads):
        if len(reads)!=2:
            if len(reads)==1:
                self.setRejectionReason(reads[0], 'unmapped mate')
            else:
                self.setRejectionReason(reads[0], 'nopair')
            return None # Only made for mate pair
        R1,R2 = reads

        self.addAlleleInfo([read for read in reads if read is not None])

        """ Valid configs:
        CATG######## R1 ########## ^ ########## R2 ##########
        ############ R2 ########## ^ ########### R1 #####CATG  reverse case
        !BWA inverts the query sequence if it maps to the negative strand!

        or R2.is_unmapped:
            if R1.is_unmapped and R2.is_unmapped:
                self.setRejectionReason(R1, 'unmapped R1;R2')
            elif R1.is_unmapped:
                self.setRejectionReason(R1, 'unmapped R1')
                self.setRejectionReason(R2, 'unmapped R1')
            else:
                self.setRejectionReason(R1, 'unmapped R2')
                self.setRejectionReason(R2, 'unmapped R2')
            return(None)
        """


        # Obtain RT hexamer:
        if R2 is not None:
            hstart, hseq = tagtools.getRandomPrimerHash(R2, onStart=True, primerLength=6)
            self.setRandomPrimer(R1,R2, hstart, hseq )

        if R1 is None or R1.is_unmapped :
            self.setRejectionReason(R1, 'unmapped R1')
            self.setRejectionReason(R2, 'unmapped R1')
            return None

        if R1.seq[:4]=='CATG' and not R1.is_reverse:
            rpos = ( R1.reference_name, R1.reference_start)
            self.addSite( [R1,R2],  strand=0, restrictionChrom=rpos[0], restrictionPos=rpos[1] )
            self.setRecognizedSequence(R1, 'CATG')
            self.setRecognizedSequence(R2, 'CATG')
            return(rpos)
        elif R1.seq[-4:]=='CATG' and R1.is_reverse:
            rpos = ( R1.reference_name, R1.reference_end-4)
            self.addSite( [R1,R2],  strand=1, restrictionChrom=rpos[0], restrictionPos=rpos[1] )
            self.setRecognizedSequence(R1, 'CATG')
            self.setRecognizedSequence(R2, 'CATG')
            return(rpos)

        # Sometimes the cycle is off
        elif R1.seq[:3]=='ATG' and not R1.is_reverse:
            rpos = ( R1.reference_name, R1.reference_start-1)
            self.addSite( [R1,R2],  strand=0, restrictionChrom=rpos[0], restrictionPos=rpos[1] )
            self.setRecognizedSequence(R1, 'ATG')
            self.setRecognizedSequence(R2, 'ATG')
            return(rpos)
        elif R1.seq[-3:]=='CAT' and R1.is_reverse: # First base was trimmed or lost
            rpos = ( R1.reference_name, R1.reference_end-3)
            self.addSite( [R1,R2],  strand=1, restrictionChrom=rpos[0], restrictionPos=rpos[1] )
            self.setRecognizedSequence(R1, 'CAT')
            self.setRecognizedSequence(R2, 'CAT')
            return(rpos)

        else:
            if R1.seq[:4]=='CATG' and R1.is_reverse:
                self.setRejectionReason(R1, 'found CATG R1 REV exp FWD')
                self.setRejectionReason(R2, 'found CATG R1 REV exp FWD')

            elif R1.seq[-4:]=='CATG' and not R1.is_reverse:
                self.setRejectionReason(R1, 'found CATG R1 FWD exp REV')
                self.setRejectionReason(R2, 'found CATG R1 FWD exp REV')
            else:
                self.setRejectionReason(R1, 'no CATG')
                self.setRejectionReason(R2, 'no CATG')
            return None

        try:
            start, end = tagtools.getPairGenomicLocations(R1,R2, R1PrimerLength=4, R2PrimerLength=6)
            self.setFragmentSize(R1, end-start)
            self.setFragmentSize(R2, end-start)
            self.setFragmentTrust(R1, start, end)
            self.setFragmentTrust(R2, start, end)

        except Exception as e:
            self.setFragmentSize(R1, 'unknown')
            self.setFragmentSize(R2, 'unknown')

        """
        if R1.seq[:4]=='CATG' and R1.reference_start<=R2.reference_start: # Site on the start of R1, R2 should map behind
            self.addSite( [R1,R2],  strand=0, restrictionChrom=R1.reference_name, restrictionPos=R1.reference_start )
            return(( R1.reference_name, R1.reference_start))

        if R1.seq[-4:]=='CATG' and R1.reference_start>=R2.reference_start: # Site on the end of R1, R2 should map before
            self.addSite( [R1,R2],  strand=1, restrictionChrom=R1.reference_name, restrictionPos=R1.reference_end-4 )
            return( (R1.reference_name, R1.reference_end-4))
        """

class ChicSeqFlagger( DigestFlagger ):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )

    def addSite(self, reads, strand, restrictionChrom, restrictionPos,is_trimmed=False ):

        sample = reads[0].get_tag(self.sampleTag)
        umi = reads[0].get_tag(self.umiTag)
        allele = None if not reads[0].has_tag(self.alleleTag) else reads[0].get_tag(self.alleleTag)
        siteInfo = tuple( [ x  for x in [strand, allele, umi] if x is not None])
        moleculeId = self.increaseAndRecordOversequencing(  sample,  restrictionChrom, restrictionPos, siteInfo=siteInfo)

        for read in reads:
            if read is None:
                continue
            self.setSiteOversequencing( read, moleculeId )
            self.setSiteCoordinate(read, restrictionPos)

            if is_trimmed:
                self.setRecognizedSequence(read,read.get_tag('lh'))
            else:
                if reads[0].is_reverse: #
                    self.setRecognizedSequence(read, reads[0].seq[-2:][::-1].translate(complement)) # the first two base, this should be A{A:80%, N:20%}, we take the complement because the reads captures the complement strand
                else:
                    self.setRecognizedSequence(read, reads[0].seq[:2]) # the last two bases
            self.setSource(read, 'CHIC')
            if allele is not None:
                self.setAllele(read, allele)

            self.setStrand(read, '+' if strand==1 else ('-' if strand==0 else '?')) #Note: strand 1 is +

    def digest(self, reads):
        if len(reads)!=2:
            return None # Only made for mate pair

        R1,R2 = reads
        if R1 is None or R2 is None:
            return None
        self.addAlleleInfo([read for read in reads if read is not None])

        """ Valid configs:
        Observed read pair:
        T######## R1 ########## ^ ########## R2 ##########
        Real molecule:
        A###########################################################

        Real molecule: and adapter:
        ADAPTER: 5' NNNNT 3'
                        A{A:80%,N:20%}NNN [CHIC MOLECULE]
                                      ^ real cut site
        """

        is_trimmed = (R1.has_tag('MX') and R1.get_tag('MX').startswith('scCHIC'))

        if R1.is_unmapped or R2.is_unmapped:
            return(None)
        try:
            start, end = tagtools.getPairGenomicLocations(R1,R2, R1PrimerLength=1 - int(is_trimmed), R2PrimerLength=6)
            self.setFragmentSize(R1, end-start)
            self.setFragmentSize(R2, end-start)
            self.setFragmentTrust(R1, start, end)
            self.setFragmentTrust(R2, start, end)

        except Exception as e:
            self.setFragmentSize(R1, 'unknown')
            self.setFragmentSize(R2, 'unknown')


        #if R1.seq[0]=='T': # Site on the start of R1, R2 should map behind
        if is_trimmed:
            # The first base of the read has been taken off and the lh tag is already set, this can be copied to RZ

            self.addSite( [R1,R2],
                strand=int(R1.is_reverse), # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                restrictionChrom=R1.reference_name,
                restrictionPos=(R1.reference_end if R1.is_reverse else R1.reference_start),
                is_trimmed=True

                  )

        else:

            self.addSite( [R1,R2],
                strand=int(R1.is_reverse), # We sequence the other strand (Starting with a T, this is an A in the molecule), the digestion thus happened on the other strand
                # On the next line we asume that the mnsase cut is one base after the ligated A, but it can be more bases upstream
                restrictionChrom=R1.reference_name,
                restrictionPos=(R1.reference_end-1 if R1.is_reverse else R1.reference_start+1),
                is_trimmed=False)
        return(( R1.reference_name, R1.reference_start))


class TagFlagger( DigestFlagger ):

    def __init__(self, tag=None,**kwargs):
        DigestFlagger.__init__(self, **kwargs )
        if tag is None or (tag!='chrom' and len(tag)!=2):
            raise ValueError(f'Invalid tag:{tag}')
        self.tag = tag


    def addSite(self, reads, strand, siteDef ):
        sample = reads[0].get_tag(self.sampleTag)
        if reads[0].has_tag(self.umiTag):
            umi = reads[0].get_tag(self.umiTag)
        else:
            umi = 'N'
        allele = None if not reads[0].has_tag(self.alleleTag) else reads[0].get_tag(self.alleleTag)
        siteInfo = tuple( [ x  for x in [strand, allele, umi,siteDef] if x is not None])


        moleculeId = self.increaseAndRecordOversequencing(
            sample,
            reads[0].reference_name,
            0, siteInfo=siteInfo)
        #print(siteDef,siteInfo,moleculeId)

        for read in reads:
            self.setSiteOversequencing( read, moleculeId )
            self.setSiteCoordinate(read, 0)
            self.setRecognizedSequence(read,f'TAG:{self.tag}')
            self.setSource(read, 'TAG')
            if allele is not None:
                self.setAllele(read, allele)

            self.setStrand(read, '+' if strand==1 else ('-' if strand==0 else '?')) #Note: strand 1 is +


    def digest(self, reads):
        self.addAlleleInfo([read for read in reads if read is not None])
        R1= reads[0]
        if len(reads)==2:
            R2=reads[1]
        else:
            R2=None
        for read in reads:
            if self.tag=='chrom':
                if R1.reference_name!=None:
                    siteDef = str(R1.reference_name)
                if R2 is not None and R2.reference_name!=None:
                    siteDef = str(R2.reference_name)
            else:
                if R1.has_tag(self.tag):
                    siteDef = R1.get_tag(self.tag)
                elif R2 is not None and R2.has_tag(self.tag):
                    siteDef = R2.get_tag(self.tag)
                else:
                    return None

            if reads[0].is_read1:
                strand = reads[0].is_reverse
            else:
                strand = not reads[0].is_reverse
            self.addSite( reads,
                strand=int(strand),
                siteDef=siteDef
            )


class MSPJIFlagger(DigestFlagger ):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )
        if self.reference is None:
            raise ValueError('The MSPJI digest flagger requires a reference genome file (FASTA)')

    def addSite(self, reads, strand, context, restrictionChrom, restrictionPos, ligationSite ):
        sample = reads[0].get_tag(self.sampleTag)
        umi = reads[0].get_tag(self.umiTag)
        allele = None if not reads[0].has_tag(self.alleleTag) else reads[0].get_tag(self.alleleTag)
        siteInfo = tuple( [ x  for x in [strand, allele, umi] if x is not None])
        moleculeId = self.increaseAndRecordOversequencing(  sample,  restrictionChrom, restrictionPos, siteInfo=siteInfo)

        for read in reads:
            self.setSiteOversequencing( read, moleculeId )
            self.setSiteCoordinate(read, restrictionPos)
            self.setRecognizedSequence(read, context)
            self.setSource(read, 'MSPJI')
            self.setLigationSite(read, ligationSite)
            if allele is not None:
                self.setAllele(read, allele)

            self.setStrand(read, strand)

    def digest(self, reads):
        if len(reads)!=2:
            return None # Only made for mate pair
        R1,R2 = reads
        if R1 is None or R2 is None:
            return None
        window = 20
        allele = self.addAlleleInfo(reads)
        drop=False
        for operation, length in R1.cigartuples:
            if operation!=0:
                drop=True
                break

        methylatedSite = None

        # CNNR NNNN NNNN | NNNN
        # CNNY NNNN NNNN  NNNNN
        if not R1.is_reverse:
            adapter = R1.seq[:4]
            seq =  self.reference.fetch( R1.reference_name, R1.reference_start-window, R1.reference_start+window+1 )
            detectedForwardFlank = seq[window-13]=='C'
            detectedForwardFlankMotif = seq[window-13:window-11]
            detectedForwardFragmentMotif = seq[window+15:window+17].translate(complement)[::-1]
            detectedForwardFragment = detectedForwardFragmentMotif[0]=='C'

            if detectedForwardFlank and not detectedForwardFragment:
                patternType = f'FWD13_{detectedForwardFlankMotif}'
                methylatedSite = R1.reference_start-13
                strand = '+'
            elif detectedForwardFragment and not detectedForwardFlank:
                patternType = f'REV16_{detectedForwardFragmentMotif}'
                methylatedSite = R1.reference_start-16
                strand = '-'

            elif detectedForwardFlank and detectedForwardFragment:
                patternType = 'FWDB'
                methylatedSite = None
            else:
                patternType = 'FWDO'
                methylatedSite = None
        else:
            adapter = R1.seq[-4:].translate(complement)[::-1]
            seq =  self.reference.fetch( R1.reference_name, R1.reference_end-window-1, R1.reference_end+window )
            detectedReverseFlank = seq[window-16]=='C'
            detectedReverseFlankMotif = seq[window-16:window-14]
            detectedReverseFragmentMotif = seq[window+12:window+14].translate(complement)[::-1]
            detectedReverseFragment = detectedReverseFragmentMotif[0]=='C'
            if detectedReverseFlank and not detectedReverseFragment:
                patternType = f'FWD16_{detectedReverseFlankMotif}'
                methylatedSite = R1.reference_start-16
                strand = '+'
            elif detectedReverseFragment and not detectedReverseFlank:
                patternType = f'REV12_{detectedReverseFragmentMotif}'
                methylatedSite = R1.reference_start-12
                strand = '-'
            elif detectedReverseFragment and detectedReverseFlank:
                methylatedSite = None
                patternType = f'REVB'
            else:
                patternType = f'REVO'
                methylatedSite = None

        #Check if we already saw this site:
        if methylatedSite is None:
            return None


        self.addSite( [R1,R2],
            strand=strand,
            restrictionChrom=R1.reference_name,
            restrictionPos=methylatedSite,
            context=patternType.split('_')[-1],
            ligationSite=adapter
         )

class ScarFlagger( DigestFlagger ):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )

    def addSite(self, reads, scarChromosome, scarPrimerStart ):

        R1_primer_length = 20
        R2_primer_length = 18

        sample = reads[0].get_tag(self.sampleTag)
        allele = None if not reads[0].has_tag(self.alleleTag) else reads[0].get_tag(self.alleleTag)

        R1, R2 = reads
        if R1 is None or R2 is None:
            for read in reads:
                if read is not None:
                    self.setRejectionReason(read,'NotPaired')

            return None
        # find all deletions:
        scarDescription = set()

        qualities = []
        for read in reads:

            firstCigarOperation, firstCigarLen = read.cigartuples[0]
            insertPos = 0
            lastNonNoneRefPos = read.reference_start if firstCigarOperation!=4 else  read.reference_start - firstCigarLen

            expandedCigar= []
            for cigarOperation, cigarLen in read.cigartuples:
                expandedCigar+=[cigarOperation]*cigarLen


            for queryPos, referencePos in read.get_aligned_pairs(matches_only=False):
                if queryPos is not None :
                    qualities.append( phred_to_prob(read.qual[queryPos]) )

                if queryPos is None and referencePos is None:
                    continue

                if referencePos is not None:
                    lastNonNoneRefPos = referencePos
                    insertPos = 0
                    if queryPos is not None: # If both the reference and query match, there is not scar information
                        continue

                if queryPos is None:
                    scarDescription.add( f'{referencePos}.D' )
                elif referencePos is None: # insert or clip:
                    operation = expandedCigar[queryPos]
                    if operation==1:
                        if lastNonNoneRefPos is None:
                            raise ValueError('Unsolvable :(')
                        queryBase = read.seq[queryPos]
                        scarDescription.add( f'{queryBase}.{insertPos+lastNonNoneRefPos}.I' )
                        insertPos+=1

        scarDescription = ','.join(sorted(list(scarDescription)))

        siteInfo = tuple( [ x  for x in [ allele, scarDescription] if x is not None])

        moleculeId = self.increaseAndRecordOversequencing(  sample,  scarChromosome, scarPrimerStart, siteInfo=siteInfo)

        # Add average base calling quality excluding primers:
        meanQ = np.mean(qualities)
        for read in reads:

            read.set_tag('SQ',1-meanQ)

            self.setSiteOversequencing( read, moleculeId )
            if len(scarDescription)==0:
                scarDescription = 'WT'
            read.set_tag('SD',scarDescription)
            self.setSiteCoordinate(read, R1.reference_start)
            self.setRecognizedSequence(read, 'SCAR')
            self.setSource(read, 'SCAR')
            if allele is not None:
                self.setAllele(read, allele)


    def digest(self, reads):
        if len(reads)!=2:
            return None # Only made for mate pair
        R1,R2 = reads

        self.addAlleleInfo(reads)
        if R1 is None or R1.is_unmapped or R2 is None or R2.is_unmapped:
            return(None)

        self.addSite( reads,  scarChromosome=R1.reference_name, scarPrimerStart = R1.reference_start )

class AlleleTagger(DigestFlagger ):

    def __init__(self, **kwargs):
        DigestFlagger.__init__(self, **kwargs )

    def digest(self, reads):
        nonNullReads = [read for read in reads if read is not None]

        self.addAlleleInfo(nonNullReads)


if __name__ == "__main__":
    # These data sources are fed to all flaggers
    flaggerArguments ={
        'reference': None if args.ref is None else pysamIterators.CachedFasta(
                                                    pysam.FastaFile(args.ref)),
        'alleleResolver': None if args.alleles is None else alleleTools.AlleleResolver(
                                args.alleles, lazyLoad=not args.loadAllelesToMem),
        'moleculeRadius': args.moleculeRadius,
        'verbose':args.verbose,
        'exon_gtf': args.exons,
        'intron_gtf': args.introns
    }

    pairedEnd = False
    flaggers = []
    qFlagger=None
    if args.ftag:
        qFlagger = QueryNameFlagger(**flaggerArguments)
        flaggers.append( qFlagger )
    if args.nla:
        flaggers.append( NlaIIIFlagger(**flaggerArguments) )
        pairedEnd=True
    if args.mspji:
        flaggers.append( MSPJIFlagger(**flaggerArguments) )
        pairedEnd=True
    if args.chic:
        flaggers.append( ChicSeqFlagger(**flaggerArguments) )
        pairedEnd=True
    if args.scar:
        flaggers.append( ScarFlagger(**flaggerArguments) )
        pairedEnd=True
    if args.tag:
        flaggers.append( TagFlagger(tag=args.tag))
    if args.atag:
        print("Running allele tagging")
        flaggers.append( AlleleTagger(**flaggerArguments))
    if args.rna:
        flaggers.append( RNA_Flagger(**flaggerArguments))


    if args.alleles is not None:
        # Check if the variant file is valid..
        if not ( args.alleles.endswith('.vcf.gz') or args.alleles.endswith('.bcf.gz') ):
            raise ValueError(f"""Please supply an indexed (bg)zipped VCF file.
            You can convert your file using: bcftools view {args.alleles} -O b -o {args.alleles}.gz;
            then index using bcftools index {args.alleles}.gz """)
        if not (os.path.exists(args.alleles+'.csi') or os.path.exists(args.alleles+'.tbi') ):
            raise ValueError(f"""Please supply an indexed (bg)zipped VCF file.
            Index using: bcftools index {args.alleles} """)

    if pairedEnd:
        print('Assuming the input is paired end')
    else:
        print('Assuming the input is single end ( Has no influence on ftag)')

    if not os.path.exists(args.o):
        try:
            os.makedirs(args.o)
        except Exception as e:
            pass

    ## Here we walk through the bamfiles and fetch

    for bamFilePath in args.bamfiles:

        bamFile = pysam.AlignmentFile(bamFilePath, "rb")
        header = bamFile.header.copy()
        outPathTemp = f'{args.o}/{os.path.basename(bamFilePath)}.unsorted'
        outPathTempWithHeader = f'{args.o}/{os.path.basename(bamFilePath)}.unsorted.with_header.bam'
        outPath = f"{args.o}/{os.path.basename(bamFilePath).replace('.bam', '')}.bam"

        print(f'Now reading {bamFilePath}, writing to {outPath}')
        if args.dedup:
            dedupOutPathTemp = f'{args.o}/{os.path.basename(bamFilePath)}.dedup.unsorted'
            dedupOutPathTempWithHeader = f'{args.o}/{os.path.basename(bamFilePath)}.dedup.with_header.bam'
            dedupOutPath = f"{args.o}/{os.path.basename(bamFilePath).replace('.bam', '')}.dedup.bam"
            dedupOutputFile = pysam.AlignmentFile(dedupOutPathTemp, "wb", header=header)
        outputFile = pysam.AlignmentFile(outPathTemp, "wb", header=header)


        bamName = os.path.basename(bamFilePath).replace('.bam','')
        #if os.path.exists(outPath+'.bai'):
        #    continue
        i=0
        if pairedEnd:
            it = pysamIterators.MatePairIterator( bamFile,
                performProperPairCheck=False, contig=args.chr)

            for i, (R1, R2) in enumerate( it ):
                for flagger in flaggers:
                    try:
                        flagger.digest( [R1,R2] )
                    except Exception as e:
                        print(e)

                        if args.fatal:
                            print(R1,R2)
                            raise e

                try:
                    if R1 is not None:
                        outputFile.write(R1)
                    if R2 is not None:
                        outputFile.write(R2)
                except Exception as e:
                    for flagger in flaggers:
                        print(flagger)
                    raise e
                if args.dedup and ((R1 is not None and R1.has_tag('RC') and R1.get_tag('RC')==1) or  (R2 is not None and R2.has_tag('RC') and R2.get_tag('RC')==1)):
                    if R1 is not None:
                        dedupOutputFile.write(R1)
                    if R2 is not None:
                        dedupOutputFile.write(R2)

                if args.head is not None and i>=args.head:
                    break
        else:
            iterator = bamFile.fetch( contig=args.chr) if args.chr is not None else bamFile
            for i, R1 in enumerate( iterator ):

                for flagger in flaggers:
                    try:
                        flagger.digest( [R1])
                    except Exception as e:
                        if args.fatal:
                            raise e
                        print(e)

                outputFile.write(R1)

                if args.head is not None and i>=args.head:
                    break
        print(f'{bamFilePath}: processed {i} reads')

        print("Finishing files (forcing close)")
        outputFile.close()
        if args.dedup:
            dedupOutputFile.close()

        if qFlagger is not None:
            readGroups={}

            for readGroup in qFlagger.assignedReadGroups:
                flowCell,lane,sampleLib = readGroup.split('.')
                library,sample = sampleLib.rsplit('_',1)
                readGroups[readGroup] = {'ID':readGroup,'LB':library, 'PL':'ILLUMINA', 'SM':sample, 'PU':readGroup}

            headerSamFilePath = outPath.replace('.bam','')+'.header.sam'
            hCopy = header.to_dict()
            hCopy['RG'] = list(readGroups.values())
            with gzip.open(outPathTemp.replace('.bam','')+'.readGroups.pickle.gz','wb') as rgzip,   pysam.AlignmentFile(headerSamFilePath,'w',header=hCopy) as headerSam:
                pickle.dump(readGroups, rgzip)

                #for readGroup, fields in readGroups.items():
                #    headerSam.write(f'@RG\tID:{readGroup}\tPL:{fields["PL"]}\tSM:{fields["SM"]}\tPU:{fields["PU"]}\tLB:{fields["LB"]}\n')
            # Clear the read groups
            qFlagger.assignedReadGroups=set()


            # Perform a reheading, sort and index
            rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {outPathTemp}; }} | samtools view -bS > {outPathTempWithHeader} ;
            rm {outPathTemp};
            samtools sort {outPathTempWithHeader} > {outPath}; samtools index {outPath};
            samtools index {outPath};
            rm {outPathTempWithHeader};
            """
            print(f"Adding read groups to header and sorting.")
            os.system(rehead_cmd)

            # Same procedure for dedup:
            if args.dedup:
                rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {dedupOutPathTemp}; }} | samtools view -bS > {dedupOutPathTempWithHeader} ;
                rm {dedupOutPathTemp};
                samtools sort {dedupOutPathTempWithHeader} > {dedupOutPath}; samtools index {dedupOutPath};
                samtools index {dedupOutPath};
                rm {dedupOutPathTempWithHeader};
                """
                print(f"Dedup file: adding read groups to header and sorting.")
                os.system(rehead_cmd)

        else:
            # we cannot assign readgroups...
            rehead_cmd = f"""
            samtools sort {outPathTemp} > {outPath}; samtools index {outPath};
            samtools index {outPath};
            rm {outPathTemp};
            """
            print(f"Adding read groups to header and sorting.")
            os.system(rehead_cmd)

            if args.dedup:
                rehead_cmd = f"""
                samtools sort {dedupOutPathTemp} > {dedupOutPath}; samtools index {dedupOutPath};
                samtools index {dedupOutPath};
                rm {dedupOutPathTemp};
                """
                print(f"Dedup file: adding read groups to header and sorting.")
                os.system(rehead_cmd)




    """
    Is:NS500414;RN:455;Fc:HYLVHBGX5;La:3;Ti:13601;CX:9882;CY:17671;Fi:N;CN:0;aa:CCGTCC;aA:CCGTCC;aI:16;LY:A3-P15-1-1;RX:ACG;RQ:GGG;BI:17;bc:ACTCGATG;BC:ACTCGATG;QT:GGKKKKKK;MX:NLAIII384C8U3;A2:TGG;AQ:E6E
    for f in $(ls *.bam | cut -f 1,2,3 -d '-' | sort | uniq | grep -v bam); do submission.py -y -time 20 -t 1 "addTagsToLennartNLABam.py $f*cell*.bam"; done

    """
