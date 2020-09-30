#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
import argparse
import collections
import functools
import gzip
import os
from singlecellmultiomics.utils import Prefetcher

def get_allele_dict():
    return collections.defaultdict(nested_set_defaultdict)
def nested_set_defaultdict():
    return collections.defaultdict( set_defaultdict )
def set_defaultdict():
    return collections.defaultdict (set)

class AlleleResolver(Prefetcher):

    def clean_vcf_name(self, vcffile):
        # Convert file name to string if it is not
        if type(vcffile) == pysam.libcbcf.VariantFile:
            vcffile = vcffile.filename.decode('ascii')
        if type(vcffile) == bytes:
            vcffile = vcffile.decode('ascii')
        return vcffile

    def __init__(self, vcffile=None,
                 chrom=None,
                 phased=True,
                 uglyMode=False,
                 lazyLoad=False,
                 select_samples=None,
                 use_cache=False,
                 ignore_conversions=None,
                 verbose=False,
                 region_start=None,
                 region_end=None


                 # When this flag is true a cache file is generated containing
                 # usable SNPs for every chromosome in gzipped format
                 ):
        """Initialise AlleleResolver

        Args:
            vcffile (str):  path of vcf file

            chrom (str):  contig/chromosome to prepare variants for

            phased(bool) : the variants in the vcf file are phased (every sample is one haplotype)

            uglyMode (bool) : the vcf file is invalid (not indexed) and has to be loaded to memory

            verbose (bool) : Print debug information

            lazyLoad (bool) : the vcf file is valid and indexed and does not need to be loaded to memory

            select_samples (list) : Use only these samples from the VCF file

            use_cache (bool) : When this flag is true a cache file is generated containing usable SNPs for every chromosome in gzipped format

            ignore_conversions(set) : conversions to ignore {(ref, alt), ..} , for example set( ('C','T'), ('G','A') )

            region_start(int) : only load variants within this range (region_start-region_end) when reading a cached variant file

            region_end(int) : only load variants within this range (region_start-region_end) when reading a cached variant file


        """
        self.args = locals().copy()
        del self.args['self']

        self.ignore_conversions = ignore_conversions
        self.phased = phased
        self.verbose = verbose
        self.locationToAllele = get_allele_dict()  # chrom -> pos-> base -> sample(s)
        self.select_samples = select_samples
        self.region_start = region_start
        self.region_end = region_end

        self.lazyLoad = lazyLoad

        if vcffile is None:
            return
        self.vcffile = self.clean_vcf_name(vcffile)
        self.use_cache = use_cache
        if use_cache:
            lazyLoad = True

        try:
            with pysam.VariantFile(vcffile) as f:
                pass
        except Exception as e:
            print(e)
            uglyMode = True
            if self.verbose:
                print(
                    'Pysam cannot read the VCF, rolling back to slower homebrew parser')

        if lazyLoad and uglyMode:
            if self.verbose:
                print(
                    'Lazy loading is not supported for non proper VCFs. Loading all variants to memory.')
            if self.select_samples is not None:
                raise NotImplementedError(
                    "Sample selection is not implemented for non proper VCF")
            lazyLoad = False

        # collections.defaultdict(set) ) #(chrom, pos)-> base -> sample(s)

        if uglyMode:
            foundIt = False  # If the target chromosome was found
            with open(vcffile, 'r') as f:

                for line in f:
                    if line[0] != '#':
                        chromosome, location, id, ref, alt = line.strip().split()
                        if chrom is None or chromosome == chrom:
                            pos = int(location)
                            foundIt = True
                            self.locationToAllele[chromosome][pos -
                                                              1][ref].add("REF")
                            self.locationToAllele[chromosome][pos -
                                                              1][alt].add("ALT")
                        elif chrom is not None and foundIt:
                            print('Stopping for speed')
                            break

            if self.verbose:
                print('Finished loading dirty vcf file')

        else:
            if not lazyLoad:
                self.fetchChromosome(vcffile, chrom)
            else:
                if self.verbose:
                    print('Lazy loading allele file')
                pass

    def addAlleleInfoOneBased(self, chromosome, location, base, alleleName):
        self.locationToAllele[chromosome][location][base].add(alleleName)

    def write_cache(self, path, chrom):
        """Write to cache file, this will make later lookups to the chromosome faster

        Args:
            path (str):  path of the cache file

            chrom (str):  contig/chromosome to write cache file for (every contig has it's own cache)
        """
        temp_path = path + '.unfinished'
        with gzip.open(temp_path, 'wt') as f:
            for position in sorted(list(self.locationToAllele[chrom].keys())):
                for base in self.locationToAllele[chrom][position]:
                    f.write(
                        f'{position}\t{base}\t{",".join(sorted(list(self.locationToAllele[chrom][position][base])))}\n')
        os.rename(temp_path, path)

    def read_cached(self, path, chrom):
        """Read cache file

        Args:
            path (str):  path of the cache file
            chrom (str):  contig/chromosome
        """
        with gzip.open(path, 'rt') as f:
            for line in f:
                position, base, samples = line.strip().split('\t', 3)
                position = int(position)
                if self.region_start is not None and position<self.region_start:
                    continue
                if self.region_end is not None and position>self.region_end:
                    break

                self.locationToAllele[chrom][position][base] = set(samples.split(','))



    def instance(self, arg_update):
        if 'self' in self.args:
            del self.args['self']
        clone = AlleleResolver(**self.args)
        return clone


    def prefetch(self, contig, start, end):

        clone = self.instance({'region_start':start, 'region_end':end})

        #print(f'Prefetching {contig}:{start}-{end}')
        try:
            self.fetchChromosome(self.vcffile, contig, True)
        except ValueError:
            # This means the chromosome is not available
            pass
        return clone

    def fetchChromosome(self, vcffile, chrom, clear=False):
        if clear:
            self.locationToAllele = get_allele_dict()  # chrom -> pos-> base -> sample(s)

        vcffile = self.clean_vcf_name(vcffile)
        # allocate:

        # Decide if this is an allele we would may be cache?
        write_cache_file_flag = False
        if self.use_cache:
            if (chrom.startswith('KN') or chrom.startswith('KZ') or chrom.startswith(
                    'chrUn') or chrom.endswith('_random') or 'ERCC' in chrom):
                write_cache_file_flag = False
            else:
                write_cache_file_flag = True

        if self.use_cache and write_cache_file_flag:
            allele_dir = f'{os.path.abspath(vcffile)}_allele_cache/'
            if not os.path.exists(allele_dir):
                os.makedirs(allele_dir)
            cache_file_name = f'{allele_dir}/{chrom}'
            if self.select_samples is not None:
                sample_list_id = '-'.join(sorted(list(self.select_samples)))
                cache_file_name = cache_file_name + '_' + sample_list_id
            cache_file_name += '.tsv.gz'
            if os.path.exists(cache_file_name):
                if self.verbose:
                    print(f"Cached file exists at {cache_file_name}")
                self.read_cached(cache_file_name, chrom)
                return
            if self.verbose:
                print(
                    f"Cache enabled, but file is not available, creating cache file at {cache_file_name}")

        self.locationToAllele[chrom][-1]['N'].add('Nop')

        unTrusted = []
        added = 0
        if self.verbose:
            print(f'Reading variants for {chrom} ', end='')
        with pysam.VariantFile(vcffile) as v:
            try:
                for rec in v.fetch(chrom, start=self.region_start, stop=self.region_end):
                    used = False
                    bad = False
                    bases_to_alleles = collections.defaultdict(
                        set)  # base -> samples

                    if self.phased:  # variants are phased, assign a random allele

                        if len(rec.samples)==0: # File without samples

                            bases_to_alleles[rec.ref]=set('r')
                            bases_to_alleles[rec.alts[0]]=set('a')
                            self.locationToAllele[rec.chrom][rec.pos - 1]=bases_to_alleles

                        else:
                            samples_assigned = set()
                            most_assigned_base = 0
                            monomorphic=False
                            for sample, sampleData in rec.samples.items():

                                if self.select_samples is not None and sample not in self.select_samples:
                                    continue
                                for base in sampleData.alleles:
                                    if base is None:
                                        # This site is monomorphic:
                                        monomorphic=True
                                        continue
                                    if len(base) == 1:
                                        bases_to_alleles[base].add(sample)
                                        used = True
                                        samples_assigned.add(sample)
                                    else:  # This location cannot be trusted:
                                        bad = True
                            # We can prune this site if all samples are associated
                            # with the same base
                            if self.select_samples is not None and used:
                                if len(samples_assigned) != len(
                                        self.select_samples):
                                    # The site is not informative
                                    bad = True
                            if monomorphic and len(bases_to_alleles)>0:
                                bad=False
                            elif len(bases_to_alleles) < 2:
                                bad = True
                                # The site is not informative
                    else:  # not phased
                        if not all(
                                len(allele) == 1 for allele in rec.alleles):  # only select SNVs
                            bad = True
                        else:
                            bad = False
                            for allele, base in zip('UVWXYZ', rec.alleles):
                                bases_to_alleles[base].add(allele)
                                used = True

                    if not bad and self.ignore_conversions is not None:  # prune conversions which are banned
                        bad = any(
                            ((rec.ref, base) in self.ignore_conversions for base in bases_to_alleles))

                    if used and not bad:
                        self.locationToAllele[rec.chrom][rec.pos -
                                                         1] = bases_to_alleles
                        added += 1
            except Exception as e:
                raise
        # for t in unTrusted:
        #    if t in self.locationToAllele:
        #        del self.locationToAllele[t[0]][t[1]]
        #del unTrusted
        if self.verbose:
            print(f'{added} variants [OK]')
        if self.use_cache and write_cache_file_flag:
            if self.verbose:
                print("writing cache file")
            try:
                self.write_cache(cache_file_name, chrom)
            except Exception as e:
                if self.verbose:
                    print(f"Exception writing cache: {e}")
                pass  # @todo

    def getAllele(self, reads):
        alleles = set()
        for read in reads:
            if read is None or read.is_unmapped:
                continue
            chrom = read.reference_name
            for readPos, refPos in read.get_aligned_pairs(matches_only=True):
                readBase = read.query_sequence[readPos]
                c = self.getAllelesAt(chrom, refPos, readBase)
                if c is not None and len(c) == 1:
                    alleles.update(c)
        return alleles

    # @functools.lru_cache(maxsize=1000) not necessary anymore... complete data is already saved in dict
    def has_location(self, chrom ,pos):
        if self.lazyLoad and chrom not in self.locationToAllele:
            try:
                self.fetchChromosome(self.vcffile, chrom, clear=True)
            except Exception as e:
                if 'invalid contig' in str(e):
                    # The contig does not exists
                    return True
                if 'fetch requires an index' in str(e):
                    raise Exception('The variant file used for allele resolving does not have an index file. Use bcftools index, or vcftools index to generate an index')
                print('ERROR, in has_location (Allele Resolver):', e)
                pass
        if chrom not in self.locationToAllele or pos not in self.locationToAllele[chrom]:
            return False
        return True


    def getAllelesAt(self, chrom, pos, base):
        if self.lazyLoad and chrom not in self.locationToAllele:
            try:
                self.fetchChromosome(self.vcffile, chrom, clear=True)
            except Exception as e:
                if 'fetch requires an index' in str(e):
                    raise Exception('The variant file used for allele resolving does not have an index file. Use bcftools index, or vcftools index to generate an index')

                print('ERROR, in getAllelesAt:', e)
                pass


        if chrom not in self.locationToAllele or pos not in self.locationToAllele[chrom]:
            return None

        if base not in self.locationToAllele[chrom][pos]:
            return None

        return self.locationToAllele[chrom][pos][base]

# FWD_-13_C REV_-16_C
# FWD_-16_C REV_+12_C
