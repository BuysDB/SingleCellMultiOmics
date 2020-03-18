#!/usr/bin/env python3

import pysam
from collections import Counter
import argparse

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Alter heterozygous SNPs in the reference vcf to the SNPs present in your data. Currently only works for 129Sv/Bl6 data.""")
    argparser.add_argument('-bamfile', '-bam', metavar='bamfile', help = "Input bam file. This is used to detect which heterozygous SNPs are present in the library.",
                           type=str, required=True)
    argparser.add_argument('-vcffile','-vcf', metavar='vcffile', help = "Input vcf file. This file will be altered based on the SNPS present in the supplied library.",
                           type=str, required=True)
    argparser.add_argument('-allele', help = "define the allele that you wish to alter. Default is 129S1_SvImJ", type = str, default = "129S1_SvImJ" )
    argparser.add_argument('-o', type=str, help="Output vcf file, containing the altered heterozygous SNPs.", required=True)

    args = argparser.parse_args()

origin_vcf = pysam.VariantFile(args.vcffile)
origin_bam = pysam.AlignmentFile(args.bamfile)

### trying to alter only the genotype part of the VCF - so the 1 and 0 that assign the ref/alt base

with pysam.VariantFile(args.o, mode = 'w', header = origin_vcf.header) as output_vcf: #change to 'wb' if output is vcf.gz

    for record in origin_vcf: #.fetch('12',114682700, 114683000): #load in only one small region to check
        record.chrom, record.pos
        
        base_obs = Counter()
        pos = record.pos - 1 
        
# calculate number of observations per base on this position, determined by the input .bam file
        for pileupcolumn in origin_bam.pileup(record.chrom, pos, pos + 1):
            if pileupcolumn.pos != pos: 
                continue
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base_obs[pileupread.alignment.query_sequence[pileupread.query_position]] +=1
                    
# we need to couple the genotype (0:ref, 1, 2 etc alternative bases) to the base itself (A,C,G,T) to make the genotype compatible with the base counter
        allele_genotype = {} 
        for number in range(len(record.samples[args.allele]['GT'])):
            allele_genotype[record.samples[args.allele]['GT'][number]] = record.samples[args.allele].alleles[number]

# select which alleles to change: we only want to alter the allele if it is not reference and count is >30%    
        new_GT = tuple(GT for GT in record.samples[args.allele]['GT'] if GT != 0 # if GT = 0, it is reference
                       and (sum(base_obs.values()) > 0 # If there are no observations, nothing should be changed
                            and float(base_obs[allele_genotype[GT]])/float(sum(base_obs.values())) > 0.3)) # only change if we have >30% observations
 
 # if the length is 0, nothing should change        
        if len(new_GT) == 0:
            new_GT = record.samples[args.allele]['GT']
 
 # if length is 1, the site should be homozygous               
        if len(new_GT) == 1:
            new_GT = tuple( (new_GT[0], new_GT[0]))

        if record.samples[args.allele].alleles == (None, None):
            continue

# to prevent crashing on monomorphic sites
        if record.alts == None:
            continue

# write altered SNPs to record and to output file      
        record.samples[args.allele]['GT'] = new_GT
        output_vcf.write(record)        
