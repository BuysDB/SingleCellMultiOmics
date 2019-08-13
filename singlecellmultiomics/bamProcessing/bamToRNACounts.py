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
import multiprocessing


matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.dpi'] = 160

import scanpy as sc

def get_gene_id_to_gene_name_conversion_table(annotation_path_exons,
                                            featureTypes=['gene_name']):
    """Create a dictionary converting a gene id to other gene features,
        such as gene_name/gene_biotype etc.

    Arguments:
        annotation_path_exons(str) : path to GTF file (can be gzipped)
        featureTypes(list) : list of features to convert to, for example ['gene_name','gene_biotype']

    Returns:
        conversion_dict(dict) : { gene_id : 'firstFeature_secondFeature'}
        """
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
            # determine the conversion name:
            if 'gene_id' in keyValues and any(
                [feat in keyValues for feat in featureTypes]):
                conversion_table[keyValues['gene_id']] = '_'.join([
                    keyValues.get(feature,'None')
                    for feature in featureTypes])

    return conversion_table

def count_transcripts(cargs):
    args,contig = cargs
    if args.alleles is not None:
        allele_resolver = alleleTools.AlleleResolver(args.alleles, lazyLoad=(not args.loadAllelesToMem))
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

    # Load features
    contig_mapping=None
    #conversion_table = get_gene_id_to_gene_name_conversion_table(args.gtfexon)
    features = singlecellmultiomics.features.FeatureContainer()
    if contig_mapping is not None:
        features.remapKeys = contig_mapping
    features.loadGTF(args.gtfexon,select_feature_type=['exon'],head=args.hf,contig=contig)
    features.loadGTF(args.gtfintron,select_feature_type=['intron'],head=args.hf,contig=contig)

    # What is used for assignment of molecules?
    if args.method=='nla':
        moleculeClass = singlecellmultiomics.molecule.AnnotatedNLAIIIMolecule
        fragmentClass = singlecellmultiomics.fragment.NLAIIIFragment
        pooling_method = 1 # all data from the same cell can be dealt with separately
        stranded = None # data is not stranded
    elif args.method=='vasa' or args.method=='cs':
        moleculeClass = singlecellmultiomics.molecule.VASA
        fragmentClass = singlecellmultiomics.fragment.SingleEndTranscript
        pooling_method = 1
        stranded=1 # data is stranded
    else:
        raise ValueError("Supply a valid method")

    #COUNT:
    exon_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount
    intron_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount
    junction_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount

    gene_set = set()
    sample_set = set()
    annotated_molecules = 0
    read_molecules=0
    if args.producebam:
        bam_path_produced = f'{args.o}/output_bam_{contig}.unsorted.bam'
        with pysam.AlignmentFile(args.alignmentfiles[0]) as alignments:
            output_bam = pysam.AlignmentFile(temp_out , "wb",header=alignments.header)

    for alignmentfile_path in args.alignmentfiles:

        i=0
        with pysam.AlignmentFile(alignmentfile_path) as alignments:
            molecule_iterator = MoleculeIterator(
                alignments=alignments,
                check_eject_every=5000,
                moleculeClass= moleculeClass,
                molecule_class_args={
                    'features':features,
                    'stranded':stranded,
                    'min_max_mapping_quality':args.minmq
                },

                fragmentClass=fragmentClass,
                fragment_class_args={
                    'umi_hamming_distance':args.umi_hamming_distance,
                    'R1_primer_length':4,
                    'R2_primer_length':6},
                perform_qflag=True, # when the reads have not been tagged yet, this flag is very much required
                pooling_method=pooling_method,
                contig=contig


            )
            for i,molecule in enumerate(molecule_iterator):
                if not molecule.is_valid():
                    continue
                molecule.annotate(args.annotmethod)
                hits = molecule.hits.keys()
                allele= None
                if allele_resolver is not None:
                    allele = molecule.get_allele(allele_resolver)
                    if len(allele)==1:
                        allele = list(allele)[0]
                    else:
                        allele = 'noAllele'


                molecule.set_intron_exon_features()

                # Obtain introns/exons/splice junction information:
                for intron in molecule.introns:
                    if allele is not None:
                        intron_counts_per_cell[molecule.sample][f'{allele}_{intron}'] += 1
                    else:
                        intron_counts_per_cell[molecule.sample][intron] += 1
                for exon in molecule.exons:
                    if allele is not None:
                        intron_counts_per_cell[molecule.sample][f'{allele}_{exon}'] += 1
                    else:
                        exon_counts_per_cell[molecule.sample][exon] += 1
                for junction in molecule.junctions:
                    if allele is not None:
                        intron_counts_per_cell[molecule.sample][f'{allele}_{junction}'] += 1
                    else:
                        junction_counts_per_cell[molecule.sample][junction]+=1

                annotated_molecules += 1


                if args.head and i>args.head:
                    print(f"-head was supplied, {i} molecules discovered, stopping")
                    break
        read_molecules+=i

    if args.producebam:
        output_bam.close()


    return (
        gene_set,
        sample_set,
        junction_counts_per_cell,
        exon_counts_per_cell,
        intron_counts_per_cell,
        annotated_molecules,
        read_molecules,
        contig

        )

if __name__=='__main__':


    argparser = argparse.ArgumentParser(
     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
     description='Create count tables from BAM file.')
    argparser.add_argument('-o',  type=str, help="output data folder", default='./rna_counts/')
    argparser.add_argument('alignmentfiles',  type=str, nargs='+')
    argparser.add_argument('-gtfexon',  type=str, required=True, help="exon GTF file containing the features to plot")
    argparser.add_argument('-gtfintron',  type=str, required=True, help="intron GTF file containing the features to plot")
    argparser.add_argument('-umi_hamming_distance',  type=int, default=1)
    argparser.add_argument('-contigmapping',  type=str, help="Use this when the GTF chromosome names do not match the ones in you bam file" )
    argparser.add_argument('-minmq',  type=int, help="Minimum molcule mapping quality", default=20 )
    argparser.add_argument('-method',  type=str, help="Data type: vasa,nla,cs", required=True )
    argparser.add_argument('-head',  type=int, help="Process this amount of molecules and export tables, also set -hf to be really fast" )
    argparser.add_argument('-hf',  type=int, help="headfeatures Process this amount features and then continue, for a quick test set this to 1000 or so." )
    argparser.add_argument('-alleles',  type=str, help="Allele file (VCF)" )
    argparser.add_argument('--loadAllelesToMem',  action='store_true',help='Load allele data completely into memory')
    argparser.add_argument('--producebam',  action='store_true',help='Produce bam file with counts tagged')
    argparser.add_argument('--ignoreMT',  action='store_true',help='Ignore mitochondria')
    argparser.add_argument('-t',  type=int, default=8, help="Amount of chromosomes processed in parallel" )

    argparser.add_argument('-annotmethod',  type=int, default=1, help="Annotation resolving method. 0: molecule consensus aligned blocks. 1: per read per aligned base" )




    #argparser.add_argument('-tagged_bam_out',  type=str, help="Output bam file" )


    args = argparser.parse_args()
    workers = multiprocessing.Pool(args.t)


    if not os.path.exists(args.o):
        os.makedirs(args.o)

    jobs = []
    contigs_todo = []
    with pysam.AlignmentFile(args.alignmentfiles[0]) as g:
        # sort by size and do big ones first.. this will be faster
        for _,chrom in sorted(list(zip(g.lengths,g.references)), reverse=True):

            if chrom.startswith('ERCC') or chrom.startswith('chrUn') or chrom.endswith('_random') or chrom.startswith('GL')  or chrom.startswith('JH'):
                continue
            if args.ignoreMT and chrom in ('mt','çhrMT','MT'):
                print("Ignoring mitochondria")
                continue
            jobs.append( (args, chrom) )
            contigs_todo.append(chrom)


    exon_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount
    intron_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount
    junction_counts_per_cell = collections.defaultdict(collections.Counter) # cell->gene->umiCount
    gene_set = set()
    sample_set = set()
    read_molecules = 0
    annotated_molecules = 0
    for  (
        result_gene_set,
        result_sample_set,
        result_junction_counts_per_cell,
        result_exon_counts_per_cell,
        result_intron_counts_per_cell,
        result_annotated_molecules,
        result_read_molecules,
        result_contig
        ) in workers.imap_unordered( count_transcripts, jobs ):
        # Update all:
        gene_set.update(result_gene_set)
        sample_set.update(result_sample_set)

        for cell, counts in result_junction_counts_per_cell.items():
            junction_counts_per_cell[cell].update(counts)
        for cell, counts in result_exon_counts_per_cell.items():
            exon_counts_per_cell[cell].update(counts)
        for cell, counts in result_intron_counts_per_cell.items():
            intron_counts_per_cell[cell].update(counts)
        read_molecules+=result_read_molecules
        annotated_molecules+=result_annotated_molecules
        # Now we finished counting
        contigs_todo = [x for x in contigs_todo if x!=result_contig]
        print(f'Finished {result_contig}, so far found {read_molecules} molecules, annotated {annotated_molecules}')
        print(f"Remaining contigs:{','.join(contigs_todo)}")

        print('writing current matrices')

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

        try:
            # Write scanpy file vanilla
            adata = sc.AnnData(
                complete_matrix
            )
            adata.var_names = gene_order
            adata.obs_names = sample_order
            adata.write(f'{args.o}/scanpy_vanilla.h5ad')

            # Write scanpy file, with introns
            adata = sc.AnnData(
                complete_matrix,
                layers={
                'spliced':  sparse_intron_matrix,
                'unspliced': sparse_exon_matrix
                #'junction' : sparse_junction_matrix
               }
            )
            adata.var_names = gene_order
            adata.obs_names = sample_order
            adata.write(f'{args.o}/scanpy_complete.h5ad')
        except Exception as e:
            print("Could not (yet?) write the scanpy files, error below")
            print(e)

    print("Writing final tables to dense csv files")
    pd.DataFrame(sparse_intron_matrix.todense(), columns=gene_order, index=sample_order).to_csv(f'{args.o}/introns.csv.gz' )
    pd.DataFrame(sparse_exon_matrix.todense(), columns=gene_order, index=sample_order).to_csv(f'{args.o}/exons.csv.gz' )
    pd.DataFrame(sparse_junction_matrix.todense(), columns=gene_order, index=sample_order).to_csv(f'{args.o}/junctions.csv.gz' )

    # Write as plaintext:
    adata.to_df().to_csv(f'{args.o}/counts.csv' )
