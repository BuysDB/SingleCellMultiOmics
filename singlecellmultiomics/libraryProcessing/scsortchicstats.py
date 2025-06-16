#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import os
import argparse
import pysam
import json
from multiprocessing import Pool
from collections import Counter, defaultdict
import numpy as np

def safe_read(read):
    if read.is_unmapped or \
       not read.is_proper_pair or \
       read.is_qcfail:
        return False
    if not read.is_read1:
        return False
    if read.mapping_quality<60:
        return False
    if read.has_tag('NM') and read.get_tag('NM')>2:
        return False
    # Check for indels
    if 'I' in read.cigarstring or 'S' in read.cigarstring or 'D' in read.cigarstring:
        return False
    return True



def count_reads(args):
    bam, contig = args
    passreads = 0
    unmapped_reads, qcfailreads, duplicate_reads, safe_reads = 0, 0, 0, 0
    with pysam.AlignmentFile(bam) as al:
        for read in al.fetch(contig):
            if read.is_supplementary or read.is_secondary:
                continue
            if read.is_unmapped:
                unmapped_reads+=1
            elif read.is_qcfail:
                qcfailreads+=1
            elif read.is_duplicate:
                duplicate_reads+=1
            else:
                passreads+=1
                if safe_read(read):
                    safe_reads+=1
    return contig, passreads, qcfailreads, duplicate_reads, unmapped_reads, safe_reads

def count_reads_sc(args):
    bam, contig = args
    passreads, qcfailreads, duplicate_reads = Counter(),Counter(),Counter()
    with pysam.AlignmentFile(bam) as al:
        for read in al.fetch(contig):
            sample = read.get_tag('SM')
            if read.is_qcfail:
                qcfailreads[sample]+=1
            elif read.is_duplicate:
                duplicate_reads[sample]+=1
            else:
                passreads[sample]+=1
    return contig, passreads, qcfailreads, duplicate_reads

def np_encoder(object):
    if isinstance(object, np.generic):
        return object.item()

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Obtain statistics for scSort-ChIC libraries')
    argparser.add_argument('bam',type=str, help='Path to tagged bam file')

    argparser.add_argument('featurebeds',type=str, nargs='*', help='Path to feature coverage bed files')

    argparser.add_argument('-o',type=str, help='Path to output json file')

    argparser.add_argument('-contaminants_species',  type=str,help='Path to contaminant specie files', nargs='*')

    argparser.add_argument('-t',type=int, help='Threadcount')
    argparser.add_argument('--sc',action= 'store_true', help='Generate also stats on a single cell basis')
    argparser.add_argument('-mito_id', help='Add MT / chrM here to flag the mitochondria as contaminant')

    args = argparser.parse_args()

    if args.contaminants_species is not None:
        print('Passed:', args.contaminants_species)
        contaminant_contigs = []
        contig_to_spec=dict()
        for path in args.contaminants_species:
            print(path)
            species_table = pd.read_csv(path, delimiter='\t', header=None,
                                        index_col=None)
            species_table.columns = ['species', 'contig', 'gi']
            species_table = species_table.set_index('contig')
            print(species_table)
            contig_to_spec.update( species_table['species'].to_dict() )

            contaminant_contigs += [
                (spec,contig) for contig,spec in contig_to_spec.items()
            ]
            print(contig_to_spec)
        if args.mito_id is not None:
            contaminant_contigs+= [ ('Mitochondrial', args.mito_id) ]
            contig_to_spec[args.mito_id] = 'Mitochondrial'
    else:
        contaminant_contigs = [('Cutibacterium','CP025935.1'),
         ('E. coli RHB09-C15','CP057942.1'),
         ('E. coli K-12','NC_000913.3'),
         ('E. coli Lambda Phage','J02459.1'),
         ('Mitochondrial',args.mito_id)
        ]
    # df_alignments = pd.DataFrame([ [(int(x) if i>0 else x) for i,x in enumerate(r.split('\t'))] for r in pysam.idxstats(args.bam).split('\n')])
    # df_alignments = df_alignments.set_index(0)

    cnts = dict() # pass read counts
    qcfailreads, duplicatereads, usablereads, unmapped_reads, n_passreads, n_safe_reads = 0, 0, 0, 0, 0, 0
    # The safe reads are pass reads which map uniquely, and in a proper pair

    with pysam.AlignmentFile(args.bam) as al:
        with Pool(args.t) as workers:
            for contig, passreads, _qcfailreads, _duplicatereads, _unmapped_reads, _safe_reads in workers.imap_unordered(
                    count_reads, ((args.bam,contig) for contig in list(al.references)+['*'])):
                # Only use uniquely mappable reads to determine contig counts:
                cnts[contig] = _safe_reads #dict(workers.imap_unordered( count_reads, ((args.bam,contig) for contig in al.references)))
                qcfailreads += _qcfailreads
                duplicatereads += _duplicatereads
                unmapped_reads += _unmapped_reads
                if contig not in contig_to_spec:
                    usablereads += passreads
                n_safe_reads+=_safe_reads
                n_passreads+= passreads
    cnts = pd.Series( cnts )
    safely_mappables = cnts.sum()

    totals = {}
    for path in args.featurebeds:
        df = pd.read_csv(path,delimiter='\t',header=None)
        df.columns=['contig','start','end','obs']
        totals[path.split('_')[-2]] = df['obs'].sum()

    scaffold_coverages_per_species = Counter()
    for spec, contig in contaminant_contigs:
        scaffold_coverages_per_species[spec] += cnts.get(contig,0)
    # Normalize:

    scaffold_coverage = pd.Series({
        name:float(100*n_reads/safely_mappables) for name, n_reads in scaffold_coverages_per_species.items()
    }).to_dict()


    d = {
        'feature_coverage': {k:float(v) for k,v in (100*pd.Series(totals)/usablereads).to_dict().items()},
        'scaffold_coverage':scaffold_coverage,
        'passreads':int(n_passreads),
        'duplicatereads': duplicatereads,
        'qcfailreads': qcfailreads,
        'totalreads': n_passreads+duplicatereads+qcfailreads+unmapped_reads,
        'usablereads':usablereads,
        'unmapped_reads':unmapped_reads
        }

    if args.sc:
        sc_cnts = defaultdict(Counter) #contig -> cell -> obs
        sc_qcfailreads, sc_duplicatereads = Counter(),Counter()
        with pysam.AlignmentFile(args.bam) as al:
            with Pool(args.t) as workers:
                for contig, passreads, _qcfailreads, _duplicatereads in workers.imap_unordered( count_reads_sc, ((args.bam,contig) for contig in al.references)):
                    sc_cnts[contig] += passreads #dict(workers.imap_unordered( count_reads, ((args.bam,contig) for contig in al.references)))
                    sc_qcfailreads += _qcfailreads
                    sc_duplicatereads += _duplicatereads

        sc_scaffold_cov = defaultdict(Counter) # contaminant_group_name -> cell -> count
        for contaminant_group_name, contig in contaminant_contigs:
            for cell, n in sc_cnts[contig].items():
                sc_scaffold_cov[contaminant_group_name][cell] += n


        d['sc_se_scaffold_cov'] = sc_scaffold_cov

        # Calculate amount of reads per cell:
        reads_per_cell = Counter()
        for contig, cell_obs in sc_cnts.items():
            reads_per_cell+=cell_obs

        d['sc_se_reads'] = {k:int(v) for k,v in reads_per_cell.items()}

    with open(args.o,'w') as o:
        json.dump(d, o, indent=2,default=np_encoder)
