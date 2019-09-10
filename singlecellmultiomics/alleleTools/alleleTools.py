#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
import argparse
import collections
import functools
import gzip
import os

class AlleleResolver:

    def clean_vcf_name(self, vcffile):
        # Convert file name to string if it is not
        if type(vcffile)==pysam.libcbcf.VariantFile:
            vcffile = vcffile.filename.decode('ascii')
        if type(vcffile)==bytes:
            vcffile = vcffile.decode('ascii')
        return vcffile

    def __init__(self, vcffile=None,
        chrom=None,
        uglyMode=False,
        lazyLoad=False,
        select_samples=None,
        use_cache = False # When this flag is true a cache file is generated containing usable SNPs for every chromosome in gzipped format
        ):



        self.select_samples = select_samples
        self.vcffile=self.clean_vcf_name(vcffile)
        self.lazyLoad = lazyLoad
        self.verbose = False
        """Initialise AlleleResolver

        Args:
            vcffile (str):  path of vcf file

            chrom (str):  contig/chromosome to prepare variants for

            uglyMode (bool) : the vcf file is invalid (not indexed) and has to be loaded to memory

            lazyLoad (bool) : the vcf file is valid and indexed and does not need to be loaded to memory

            select_samples (list) : Use only these samples from the VCF file

        """


        self.use_cache = use_cache
        try:
            with  pysam.VariantFile(vcffile) as f:
                pass
        except Exception as e:
            print(e)
            uglyMode=True
            if self.verbose:
                print('Pysam cannot read the VCF, rolling back to slower homebrew parser')

        if lazyLoad and uglyMode:
            if self.verbose:
                print('Lazy loading is not supported for non proper VCFs. Loading all variants to memory.')
            if self.select_samples is not None:
                raise NotImplementedError("Sample selection is not implemented for non proper VCF")
            lazyLoad = False
        #self.locationToAllele = collections.defaultdict( lambda: collections.defaultdict(set) ) #(chrom, pos)-> base -> sample(s)
        self.locationToAllele = collections.defaultdict(  lambda: collections.defaultdict(  lambda: collections.defaultdict(set) )) #chrom -> pos-> base -> sample(s)
        if vcffile is None:
            return
        if uglyMode:
            foundIt = False # If the target chromosome was found
            with open(vcffile,'r') as f:

                for line in f:
                    if line[0]!='#':
                        chromosome, location, id, ref, alt = line.strip().split()
                        if chrom is None or chromosome==chrom:
                            pos = int(location)
                            foundIt=True
                            self.locationToAllele[chromosome][pos-1][ref].add("REF")
                            self.locationToAllele[chromosome][pos-1][alt].add("ALT")
                        elif chrom is not None and foundIt:
                            print('Stopping for speed')
                            break

            if self.verbose:
                print('Finished loading dirty vcf file')

        else:
            if not lazyLoad:
                self.fetchChromosome(vcffile,chrom)
            else:
                if self.verbose:
                    print('Lazy loading allele file')
                pass

    def addAlleleInfoOneBased( self, chromosome, location, base, alleleName ):
        self.locationToAllele[chromosome][location][base].add(alleleName)

    def write_cache(self, path, chrom):
        """Write to cache file, this will make later lookups to the chromosome faster

        Args:
            path (str):  path of the cache file

            chrom (str):  contig/chromosome to write cache file for (every contig has it's own cache)
        """
        temp_path = path+'.unfinished'
        with gzip.open(temp_path,'wt') as f:
            for position in sorted(list(self.locationToAllele[chrom].keys())):
                for base in self.locationToAllele[chrom][position]:
                    f.write(f'{position}\t{base}\t{",".join(sorted(list(self.locationToAllele[chrom][position][base])))}\n')
        os.rename(temp_path,path)

    def read_cached(self,path,chrom):
        """Read cache file

        Args:
            path (str):  path of the cache file
            chrom (str):  contig/chromosome
        """
        with gzip.open(path,'rt') as f:
            for line in f:
                position, base, samples = line.strip().split('\t',3)
                self.locationToAllele[chrom][int(position)][base] = set(samples.split(','))

    def fetchChromosome(self, vcffile, chrom, clear=False):
        if clear:
            self.locationToAllele = collections.defaultdict(
                lambda: collections.defaultdict(
                lambda: collections.defaultdict(set) )) #chrom -> pos-> base -> sample(s)

        vcffile = self.clean_vcf_name(vcffile)
        #allocate:

        # Decide if this is an allele we would may be cache?
        write_cache_file_flag=False
        if self.use_cache:
            if (chrom.startswith('KN') or chrom.startswith('KZ') or chrom.startswith('chrUn') or chrom.endswith('_random') or 'ERCC' in chrom):
                write_cache_file_flag=False
            else:
                write_cache_file_flag=True

        if self.use_cache and write_cache_file_flag:
            allele_dir = f'{os.path.abspath(vcffile)}_allele_cache/'
            if not os.path.exists(allele_dir):
                os.makedirs(allele_dir)
            cache_file_name = f'{allele_dir}/{chrom}'
            if self.select_samples is not None:
                sample_list_id = '-'.join(sorted(list(self.select_samples)))
                cache_file_name = cache_file_name+'_'+sample_list_id
            cache_file_name+='.tsv.gz'
            if os.path.exists(cache_file_name):
                if self.verbose:
                    print(f"Cached file exists at {cache_file_name}")
                self.read_cached(cache_file_name,chrom)
                return
            if self.verbose:
                print(f"Cache enabled, but file is not available, creating cache file at {cache_file_name}")

        self.locationToAllele[chrom][-1]['N'].add('Nop')

        unTrusted = []
        added = 0
        if self.verbose:
            print(f'Reading variants for {chrom} ', end='')
        with pysam.VariantFile(vcffile) as v:
            for rec in v.fetch(chrom):
                used = False
                bad = False
                bases_to_alleles = collections.defaultdict(set) # base -> samples
                samples_assigned = set()
                most_assigned_base = 0
                for sample, sampleData in rec.samples.items():
                    if self.select_samples is not None and not sample in self.select_samples:
                        continue
                    for base in sampleData.alleles:
                        if base is None:
                            unTrusted.append( (rec.chrom, rec.pos ) )
                            continue
                        if len(base)==1:
                            bases_to_alleles[base].add(sample)
                            used=True
                            samples_assigned.add(sample)
                        else: # This location cannot be trusted:
                            bad = True
                # We can prune this site if all samples are associated with the same base
                if self.select_samples is not None and used:
                    if len(samples_assigned)!=len(self.select_samples):
                        # The site is not informative
                        bad=True
                if len(bases_to_alleles)<2:
                    bad=True
                    # The site is not informative

                if used and not bad:
                    self.locationToAllele[rec.chrom][ rec.pos-1] = bases_to_alleles
                    added+=1



        #for t in unTrusted:
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
                pass # @todo

    def getAllele( self, reads ):
        alleles = set()
        for read in reads:
            if read is None or read.is_unmapped:
                continue
            chrom = read.reference_name
            for  readPos, refPos in read.get_aligned_pairs(matches_only=True):
                readBase = read.query_sequence[readPos]
                c = self.getAllelesAt(chrom, refPos, readBase)
                if c is not None and len(c)==1:
                    alleles.update(c)
        return alleles


    #@functools.lru_cache(maxsize=1000) not necessary anymore... complete data is already saved in dict
    def getAllelesAt(self, chrom, pos, base):
        if self.lazyLoad and not chrom in self.locationToAllele:
            try:
                self.fetchChromosome(self.vcffile, chrom, clear=True)
            except Exception as e:
                print(e)
                pass

        if not chrom in self.locationToAllele or not pos in self.locationToAllele[chrom] :
            return None

        if not base in self.locationToAllele[chrom][pos]:
            return None

        return self.locationToAllele[chrom][pos][base]

# FWD_-13_C REV_-16_C
# FWD_-16_C REV_+12_C
