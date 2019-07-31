#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pysam
import collections
import argparse
import gzip
import pickle
import matplotlib
import numpy as np
import singlecellmultiomics
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import singlecellmultiomics.features
import pysamiterators.iterators
import pysam
import pandas as pd
import scipy.sparse
import gzip
from singlecellmultiomics.molecule import MoleculeIterator
from singlecellmultiomics.alleleTools import alleleTools

import scanpy as sc

def get_gene_id_to_gene_name_conversion_table(annotation_path_exons):
    conversion_table = {}
    with (gzip.open(annotation_path_exons,'rt') if annotation_path_exons.endswith('.gz') else open(annotation_path_exons,'r')) as t:
        for i,line in enumerate(t):
            parts = line.rstrip().split(None,8)
            keyValues = {}
            for part in parts[-1].split(';'):
                kv = part.strip().split()
                if len(kv)==2:
                    key = kv[0]
                    value = kv[1].replace('"', '')
                    keyValues[key] =  value
            if 'gene_id' in keyValues and 'gene_name' in keyValues:
                conversion_table[keyValues['gene_id']] = keyValues['gene_name']
    return conversion_table


if __name__=='__main__':

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['figure.dpi'] = 160

    argparser = argparse.ArgumentParser(
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     description='Create count tables from BAM file.')
    argparser.add_argument('-o',  type=str, help="output data folder", default='./rna_counts/')
    argparser.add_argument('alignmentfiles',  type=str, nargs='*')
    argparser.add_argument('-gtfexon',  type=str, required=True, help="exon GTF file containing the features to plot")
    argparser.add_argument('-gtfintron',  type=str, required=True, help="intron GTF file containing the features to plot")
    argparser.add_argument('-umi_hamming_distance',  type=int, default=1)
    argparser.add_argument('-contigmapping',  type=str, help="Use this when the GTF chromosome names do not match the ones in you bam file" )
    argparser.add_argument('-method',  type=str, help="Data type: either vasa or nla" )
    argparser.add_argument('-head',  type=int, help="Process this amount of molecules and export tables" )
    argparser.add_argument('-alleles',  type=str, help="Allele file (VCF)" )

    args = argparser.parse_args()

    if args.alleles is not None:
        allele_resolver = alleleTools.AlleleResolver(args.alleles)
    else:
        allele_resolver = None

    contig_mapping=None

    if args.contigmapping=='danio':
        contig_mapping ={
            '1':'CM002885.2',
            '2':'CM002886.2',
            '3':'CM002887.2',
            '4':'CM002888.2',
            '5':'CM002889.2',

            '6':'CM002890.2',
            '7':'CM002891.2',
            '8':'CM002892.2',
            '9':'CM002893.2',
            '10':'CM002894.2',
            '11':'CM002895.2',
            '12':'CM002896.2',
            '13':'CM002897.2',
            '14':'CM002898.2',
            '15':'CM002899.2',

            '16':'CM002900.2',
            '17':'CM002901.2',
            '18':'CM002902.2',
            '19':'CM002903.2',
            '20':'CM002904.2',
            '21':'CM002905.2',
            '22':'CM002906.2',
            '23':'CM002907.2',
            '24':'CM002908.2',
            '25':'CM002909.2',
        }


    if not os.path.exists(args.o):
        os.makedirs(args.o)

    # Load features
    contig_mapping=None
    conversion_table = get_gene_id_to_gene_name_conversion_table(args.gtfexon)
    features = singlecellmultiomics.features.FeatureContainer()
    if contig_mapping is not None:
        features.remapKeys = contig_mapping
    features.loadGTF(args.gtfexon,select_feature_type=['exon'])
    features.loadGTF(args.gtfintron,select_feature_type=['intron'])

    # What is used for assignment of molecules?
    if args.method=='nla':
        moleculeClass = singlecellmultiomics.molecule.AnnotatedNLAIIIMolecule
        fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment
        pooling_method = 1
    elif args.method=='vasa':
        moleculeClass = singlecellmultiomics.molecule.VASA
        fragmentClass = singlecellmultiomics.fragment.SingleEndTranscript
        pooling_method = 1

    #COUNT:
    exon_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount
    intron_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount
    junction_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount

    gene_set = set()
    sample_set = set()
    for alignmentfile_path in args.alignmentfiles:
        with pysam.AlignmentFile(alignmentfile_path) as alignments:
            molecule_iterator = MoleculeIterator(
                alignments=alignments,
                check_eject_every=5000,
                moleculeClass= moleculeClass,
                molecule_class_args={
                    'features':features
                },

                fragmentClass=fragmentClass,
                fragment_class_args={
                    'umi_hamming_distance':args.umi_hamming_distance,
                    'R1_primer_length':4,
                    'R2_primer_length':6},
                perform_qflag=True, # when the reads have not been tagged yet, this flag is very much required
                pooling_method=pooling_method

            )
            for i,molecule in enumerate(molecule_iterator):
                molecule.annotate()
                hits = molecule.hits.keys()
                allele= None
                if allele_resolver is not None:
                    allele = molecule.get_allele(allele_resolver)
                    if len(allele)==1:
                        allele = list(allele)[0]
                    else:
                        allele = 'noAllele'
                f_hits = collections.defaultdict(collections.Counter)
                for hit in hits:
                    if hit.startswith('type:exon'):

                        gene = hit.split(',')[-1].replace('gene_id:','')
                        if allele is not None:
                            gene = f'{allele}_{gene}'
                        f_hits[gene]['exon']+=1
                    elif hit.startswith('type:intron'):
                        gene = hit.split(',')[-1].replace('gene_id:','')
                        if allele is not None:
                            gene = f'{allele}_{gene}'
                        f_hits[gene]['intron']+=1

                for gene, intron_exon_hits in f_hits.items():

                    spliced=True
                    if 'intron' in intron_exon_hits:
                        spliced=False
                        if 'exon' in intron_exon_hits:
                            junction_counts_per_cell[molecule.sample][gene]+=1

                    if not spliced:
                        exon_counts_per_cell[molecule.sample][gene] += 1
                    else:
                        intron_counts_per_cell[molecule.sample][gene] += 1

                    gene_set.add(gene)
                    sample_set.add(molecule.sample)

                if args.head and args.head>=i:
                    break
    # Now we finished counting

    # freeze order of samples and genes:
    sample_order = sorted(list(sample_set))
    gene_order = sorted(list(gene_set))

    # Construct the sparse matrices:
    sparse_intron_matrix = scipy.sparse.dok_matrix((len(sample_set),len(gene_set)),dtype=np.int64)
    #sparse_intron_matrix.setdefault(0)
    sparse_exon_matrix = scipy.sparse.dok_matrix((len(sample_set),len(gene_set)),dtype=np.int64)
    #sparse_exon_matrix.setdefault(0)
    sparse_junction_matrix = scipy.sparse.dok_matrix((len(sample_set),len(gene_set)),dtype=np.int64)

    for sample_idx,sample in enumerate(sample_order):
        if sample in exon_counts_per_cell:
            for gene, counts in exon_counts_per_cell[sample].items():
                gene_idx = gene_order.index(gene)
                sparse_exon_matrix[sample_idx, gene_idx] = counts
        if sample in intron_counts_per_cell:
            for gene, counts in intron_counts_per_cell[sample].items():
                gene_idx = gene_order.index(gene)
                sparse_intron_matrix[sample_idx, gene_idx] = counts
        if sample in junction_counts_per_cell:
            for gene, counts in junction_counts_per_cell[sample].items():
                gene_idx = gene_order.index(gene)
                sparse_junction_matrix[sample_idx, gene_idx] = counts


    # Write matrices to disk
    sparse_intron_matrix = sparse_intron_matrix.tocsc()
    sparse_exon_matrix = sparse_exon_matrix.tocsc()
    sparse_junction_matrix = sparse_junction_matrix.tocsc()
    complete_matrix = sparse_intron_matrix + sparse_exon_matrix

    scipy.sparse.save_npz(f'{args.o}/sparse_complete_matrix.npz', complete_matrix)
    scipy.sparse.save_npz(f'{args.o}/sparse_intron_matrix.npz',sparse_intron_matrix)
    scipy.sparse.save_npz(f'{args.o}/sparse_exon_matrix.npz',sparse_exon_matrix)
    scipy.sparse.save_npz(f'{args.o}/sparse_junction_matrix.npz',sparse_junction_matrix)

    # Write scanpy file
    adata = sc.AnnData(
        complete_matrix,
        layers={
        'spliced':  sparse_intron_matrix,
        'unspliced': sparse_exon_matrix,
        'junction' : sparse_junction_matrix
       }
    )
    adata.var_names = gene_order
    adata.obs_names = sample_order
    adata.write(f'{args.o}/scanpy.h5ad')
