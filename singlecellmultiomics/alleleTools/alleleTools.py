#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
import argparse
import collections
import functools

class AlleleResolver:

    def __init__(self, vcffile=None, chrom=None, uglyMode=False, lazyLoad=False):
        self.vcffile=vcffile
        self.lazyLoad = lazyLoad

        try:
            with  pysam.VariantFile(vcffile) as f:
                pass
        except Exception as e:
            print(e)
            uglyMode=True
            print('Pysam cannot read the VCF, rolling back to slower homebrew parser')

        if lazyLoad and uglyMode:
            print('Lazy loading is not supported for non proper VCFs. Loading all variants to memory.')
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
            print('Finished loading dirty vcf file')

        else:
            if not lazyLoad:
                self.fetchChromosome(vcffile,chrom)
            else:
                print('Lazy loading allele file')

    def addAlleleInfoOneBased( self, chromosome, location, base, alleleName ):
        self.locationToAllele[chromosome][location][base].add(alleleName)

    def fetchChromosome(self, vcffile, chrom, clear=False):
        if clear:
            self.locationToAllele = collections.defaultdict(  lambda: collections.defaultdict(  lambda: collections.defaultdict(set) )) #chrom -> pos-> base -> sample(s)
        #allocate:
        self.locationToAllele[chrom][-1]['N'].add('Nop')

        unTrusted = []
        added = 0
        print(f'Reading variants for {chrom} ', end='')
        with pysam.VariantFile(vcffile) as v:
            for rec in v.fetch(chrom):
                for sample, sampleData in rec.samples.items():
                    for base in sampleData.alleles:
                        if base is None:
                            unTrusted.append( (rec.chrom, rec.pos ) )
                            continue
                        if len(base)==1:
                            self.locationToAllele[rec.chrom][ rec.pos-1][base].add(sample)
                            added+=1
                        else: # This location cannot be trusted:
                            unTrusted.append( (rec.chrom, rec.pos ) )

        for t in unTrusted:
            if t in self.locationToAllele:
                del self.locationToAllele[t[0]][t[1]]
        del unTrusted
        print(f'{added} variants [OK]')

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
                pass

        if not chrom in self.locationToAllele or not pos in self.locationToAllele[chrom] :
            return None

        if not base in self.locationToAllele[chrom][pos]:
            return None

        return self.locationToAllele[chrom][pos][base]

# FWD_-13_C REV_-16_C
# FWD_-16_C REV_+12_C
