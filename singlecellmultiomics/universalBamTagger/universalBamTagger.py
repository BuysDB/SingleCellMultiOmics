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

#Modules:
from singlecellmultiomics.universalBamTagger.rna import RNA_Flagger
from singlecellmultiomics.universalBamTagger.digest import DigestFlagger
from singlecellmultiomics.universalBamTagger.nlaIII import NlaIIIFlagger
from singlecellmultiomics.universalBamTagger.scchic import ChicSeqFlagger
from singlecellmultiomics.universalBamTagger.mspjI import MSPJIFlagger
from singlecellmultiomics.universalBamTagger.scar import ScarFlagger
from singlecellmultiomics.universalBamTagger.tag import TagFlagger
from singlecellmultiomics.universalBamTagger.taps import TAPSFlagger

from singlecellmultiomics.utils.sequtils import hamming_distance,phred_to_prob

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
    tagAlgs.add_argument('--taps', action='store_true', help='Add TAPS based methylation state ')
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

    if not args.mspji and not args.nla and not args.chic and not args.ftag  and not args.rna and args.tag is None and args.atag is None and args.taps is None:
        raise ValueError('Please supply any or a combination of --ftag --nla --chic --mspji --taps')

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


def molecule_to_random_primer_dict(molecule, primer_length=6, primer_read=2): #1: read1 2: read2
    rp = collections.defaultdict(list)
    for fragment in molecule:
        if fragment[primer_read-1] is not None:
            hstart, hseq = tagtools.getRandomPrimerHash(fragment[primer_read-1], onStart=True, primerLength=6)
            rp[hstart, hseq].append(fragment)
    return rp
    #for k,p in rp.items():
    #    yield p

"""
Iterate over molecules in a bam file

Parameters
----------
alignmentfile : pysam.AlignmentFile
    file to read the molecules from

look_around_radius : int
    buffer to accumulate molecules in. All fragments belonging to one molecule should fit this radius

umi_hamming_distance : int
    Edit distance on UMI, 0: only exact match, 1: single base distance

sample_select : iterable
    Iterable of samples to only select molecules from

Yields
----------
list of molecules : list [ pysam.AlignedSegment ]
[ (R1,R2), (R1,R2) ... ]

"""
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
        if assignment_a is None or assignment_b is None:
            return False

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

    def get_fragment_chromosome(self, fragment):
        for read in fragment:
            if read is not None:
                return read.reference_name

    def __iter__(self):
        for fragment in pysamIterators.MatePairIterator( self.alignmentfile, **self.pysam_kwargs, performProperPairCheck=False ):


            if fragment[0] is None or fragment[0].reference_name is None:
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
            if self.get_fragment_chromosome(fragment)!=self.current_chromosome and self.current_chromosome is not None:
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

            self.current_chromosome = self.get_fragment_chromosome(fragment)
            self.current_position = self.localisation_function(fragment) #[0].reference_end

        # yield everything which was not yielded yet
        for molecule in self._yield_all_in_current_cache():
            yield molecule


"""
Iterate over transcripts in a bam file

Parameters
----------
alignmentfile : pysam.AlignmentFile
    file to read the molecules from

look_around_radius : int
    buffer to accumulate molecules in. All fragments belonging to one molecule should fit this radius

informative_read : int
    which read is used to define the mapping coordinate of the fragment, 1 or 2.

umi_hamming_distance : int
    Edit distance on UMI, 0: only exact match, 1: single base distance, 2 base distance ...

assignment_radius : int
    tolerance on fragment starting coordinates

sample_select : iterable
    Iterable of samples to only select molecules from

Yields
----------
list of molecules : list [ pysam.AlignedSegment ]
[ (R1,R2), (R1,R2) ... ]
"""
class TranscriptIterator(MoleculeIterator):
    def __init__(self, look_around_radius=100, informative_read=2, assignment_radius=10,**kwargs):
        MoleculeIterator.__init__(self,look_around_radius=look_around_radius,**kwargs)
        self.informative_read = informative_read
        self.assignment_radius = assignment_radius

    def assignment_function(self, fragment):
        return fragment[self.informative_read-1].get_tag('SM'),fragment[self.informative_read-1].get_tag('RX'),fragment[self.informative_read-1].is_reverse

    def localisation_function(self,fragment):
        return int((fragment[self.informative_read-1].reference_start)/self.assignment_radius)*self.assignment_radius



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
    if args.taps:
        flaggers.append( TAPSFlagger(**flaggerArguments))
        pairedEnd=True


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
            rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {outPathTemp}; }} | samtools view -b > {outPathTempWithHeader} ;
            rm {outPathTemp};
            samtools sort {outPathTempWithHeader} > {outPath}; samtools index {outPath};
            rm {outPathTempWithHeader};
            """
            print(f"Adding read groups to header and sorting.")
            os.system(rehead_cmd)

            # Same procedure for dedup:
            if args.dedup:
                rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {dedupOutPathTemp}; }} | samtools view -b > {dedupOutPathTempWithHeader} ;
                rm {dedupOutPathTemp};
                samtools sort {dedupOutPathTempWithHeader} > {dedupOutPath}; samtools index {dedupOutPath};
                rm {dedupOutPathTempWithHeader};
                """
                print(f"Dedup file: adding read groups to header and sorting.")
                os.system(rehead_cmd)

        else:
            # we cannot assign readgroups...
            rehead_cmd = f"""
            samtools sort {outPathTemp} > {outPath}; samtools index {outPath};
            rm {outPathTemp};
            """
            print(f"Adding read groups to header and sorting.")
            os.system(rehead_cmd)

            if args.dedup:
                rehead_cmd = f"""
                samtools sort {dedupOutPathTemp} > {dedupOutPath}; samtools index {dedupOutPath};
                rm {dedupOutPathTemp};
                """
                print(f"Dedup file: adding read groups to header and sorting.")
                os.system(rehead_cmd)




    """
    Is:NS500414;RN:455;Fc:HYLVHBGX5;La:3;Ti:13601;CX:9882;CY:17671;Fi:N;CN:0;aa:CCGTCC;aA:CCGTCC;aI:16;LY:A3-P15-1-1;RX:ACG;RQ:GGG;BI:17;bc:ACTCGATG;BC:ACTCGATG;QT:GGKKKKKK;MX:NLAIII384C8U3;A2:TGG;AQ:E6E
    for f in $(ls *.bam | cut -f 1,2,3 -d '-' | sort | uniq | grep -v bam); do submission.py -y -time 20 -t 1 "addTagsToLennartNLABam.py $f*cell*.bam"; done

    """
