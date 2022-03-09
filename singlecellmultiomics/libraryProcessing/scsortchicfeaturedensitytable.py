#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from singlecellmultiomics.features import FeatureContainer
from multiprocessing import Pool
from singlecellmultiomics.utils import pool_wrapper
from more_itertools import chunked
import pandas as pd
import numpy as np
import os
import argparse
import pysam
import gzip


def coord_to_index_clip(coordinate,centroid,radius,mirror,bp_per_bin):

    dist = centroid - coordinate
    if not mirror:
        dist = -dist
    index = dist

    index = int(index/bp_per_bin)
    index+=int(radius/bp_per_bin)

    return np.clip(index,0,radius*2)

def sample_locations(targets, bam_path, radius, bp_per_bin, extra_sample_rad=500):

    obs = np.zeros( (len(targets), int((radius*2+1)/bp_per_bin) ))
    with  pysam.AlignmentFile(bam_path,threads=8) as alns:

        for i,(id, (f_contig,f_start,f_end,strand)) in enumerate(targets.items()):

            if strand=='+':
                centroid= f_start
                mirror = False
            elif strand=='-':
                centroid = f_end
                mirror = True

            fetch_coords = [centroid-radius-extra_sample_rad,centroid+radius+extra_sample_rad]
            fetch_coords = np.clip(fetch_coords,0, alns.get_reference_length(f_contig))
            for read in alns.fetch(f_contig, *fetch_coords):
                if not read.has_tag('DS') or read.is_duplicate:
                    continue

        #         # Based on DS:
        #         index = coord_to_index(read.get_tag('DS'),centroid,radius,mirror)
        #         if index is not None:
        #             obs[0,index]+=1

                # Based on coordinates
                index_start= coord_to_index_clip(read.reference_start,centroid,radius,mirror,bp_per_bin)
                index_end= coord_to_index_clip(read.reference_end,centroid,radius,mirror,bp_per_bin)
                if index_start!=index_end:
                    index_start,index_end = min(index_start,index_end), max(index_start,index_end)
                    obs[i,index_start:index_end]+=1

    index = pd.MultiIndex.from_tuples( list(targets.keys()) )
    df  = pd.DataFrame(obs,index=index)
    return df





def coord_to_index(coordinate,centroid,radius,mirror):

    dist = centroid - coordinate
    if not mirror:
        dist = -dist
    index = dist
    index+=radius
    if index>=0 and dist<radius:
        return index



def get_density(bam_path, radius, bp_per_bin, targets, n_threads):
    with Pool(n_threads) as workers:



        dfs = list(workers.imap(pool_wrapper,
                                  ((sample_locations,
                                       {
                                           'targets': { k:targets[k] for k in keys },
                                           'bam_path':bam_path,
                                           'radius':radius,
                                           'bp_per_bin':bp_per_bin,
                                           'extra_sample_rad':1000
                                       }
                                        )
                                       for keys in chunked(targets.keys(),3000)
                                  )))
        df = pd.concat(dfs,axis=0)
        df.columns = (df.columns*bp_per_bin) - radius
        return df


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Obtain count table over gene bodies for scSort-ChIC libraries')
    argparser.add_argument('bam',type=str, help='Path to tagged bam file')

    argparser.add_argument('featuregtf',type=str, help='Path to feature gtf file')

    argparser.add_argument('-o',type=str, help='Path to output csv.gz file')

    argparser.add_argument('-radius',type=int, help='Radius in bp around feature',default=100_000)

    argparser.add_argument('-bp_per_bin',type=int, help='bp per bin',default=250)

    argparser.add_argument('-t',type=int, help='threads')

    argparser.add_argument('-feature',type=str, default='gene', help='Selected feature type')

    #argparser.add_argument('-feature',type=int, help='feature to be considered for counting')

    args = argparser.parse_args()

    feature = args.feature
    features = FeatureContainer()
    features.loadGTF(args.featuregtf, select_feature_type=[feature],store_all=True)

    targets = {}
    prev=None
    metas = {}
    got_genes = set() # Do not count genes more than once.
    for i,(f_contig,f_start,f_end,id,strand,meta) in enumerate(features):
        if prev==(f_contig,f_start,f_end,strand):
            continue
        if id in got_genes:
            continue
        prev=(f_contig,f_start,f_end,strand)

        meta = dict(meta)

        gene_name = meta.get('gene_name',id)
        if gene_name in got_genes:
            continue
        got_genes.add(id)
        got_genes.add(gene_name)

        targets[(strand,id)] = prev
        meta['length'] = f_end - f_start
        metas[(strand,id)] = meta


    d = get_density(args.bam, args.radius, args.bp_per_bin, targets, args.t)
    d.index = [ metas[idx].get('gene_name',idx) for idx in d.index]
    d.to_csv(args.o)
