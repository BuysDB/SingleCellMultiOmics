#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pysam
import pandas as pd
import numpy as np
from multiprocessing import Pool, Manager
from collections import defaultdict
from singlecellmultiomics.bamProcessing import get_reference_path_from_bam
from singlecellmultiomics.molecule import MoleculeIterator,TAPS

def get_bulk_vector(args):
    obj, samples, location = args
    return obj.get_bulk_column(samples, location)

class MethylationCountMatrix:

    def __init__(self, counts: dict = None, threads=None):
        # Sample->(contig,bin_start,bin_end)-> [methylated_counts, unmethylated]
        self.counts = {} if counts is None else counts

        # { (contig, bin_start, bin_end), (contig, bin_start, bin_end) .. }
        #or
        # { (contig, bin_start, bin_end,strand), (contig, bin_start, bin_end, strand) .. }
        self.sites = set()
        self.threads = threads

    def __getitem__(self, key: tuple):
        sample, location = key
        if not sample in self.counts:
            self.counts[sample] = {}
        if not location in self.counts[sample]:
            self.sites.add(location)
            self.counts[sample][location] = [0, 0]
        return self.counts[sample][location]

    def get_without_init(self, key: tuple):
        # Obtain a key without setting it
        # sample, location = key
        try:
            return self.counts[key[0]][key[1]]
        except KeyError:
            return (0,0)

    def __setitem__(self, key: tuple, value: list):
        sample, location = key
        if not sample in self.counts:
            self.counts[sample] = {}
        self.counts[sample][location] = value

    def update(self, other):
        # This does not work for regions with overlap! Those will be overwritten
        for sample, counts in other.counts.items():
            if sample not in self.counts:
                self.counts[sample] = {}
            self.counts[sample].update(counts)
        self.sites.update(other.sites)

    def get_sample_list(self):
        return sorted(list(self.counts.keys()))

    def __repr__(self):
        return f'Methylation call matrix containing {len(self.counts)} samples and {len(self.sites)} locations'

    def prune(self, min_samples: int = 0, min_variance: float = None):
        if len(self.sites)==0 or len(self.counts) == 0 or min_samples == 0 and min_variance is None:
            return

        for location, row in self.get_bulk_frame(use_multi=False).iterrows():
            if row.n_samples < min_samples:
                self.delete_location(location)
            elif min_variance is not None and (np.isnan(row.variance) or row.variance < min_variance):
                self.delete_location(location)

    def delete_location(self, location):

        drop_samples = []
        for sample in self.counts:
            if location in self.counts[sample]:
                del self.counts[sample][location]
                if len(self.counts[sample]) == 0:
                    drop_samples.append(sample)
        self.sites.remove(location)

        # Remove samples without any data left:
        for d in drop_samples:
            del self.counts[d]

    def get_sample_distance_matrix(self):
        self.check_integrity()
        def distance(row, matrix):
            # Amount of differences / total comparisons
            return np.nansum(np.abs((matrix - row)), axis=1) / (np.isfinite(matrix - row).sum(axis=1))

        def get_dmat(df):
            dmat = np.apply_along_axis(distance, 1, df.values, matrix=df.values)
            return pd.DataFrame(dmat, columns=df.index, index=df.index)

        with np.errstate(divide='ignore', invalid='ignore'):
            dmat = get_dmat(self.get_frame('beta'))

            while dmat.isna().sum().sum() > 0:
                sample = dmat.isna().sum().idxmax()
                dmat.drop(sample, 0, inplace=True)
                dmat.drop(sample, 1, inplace=True)

        return dmat




    def get_frame(self, dtype: str):
        """
        Get pandas dataframe containing the selected column

        Args:
            dtype: either 'methylated', 'unmethylated' or 'beta'

        Returns:
            df(pd.DataFrame) : Dataframe containing the selected column, rows are samples, columns are locations

        """
        self.check_integrity()
        # Fix columns
        columns = list(sorted(self.sites))
        # Create column to index mapping:
        column_to_index = {c: i for i, c in enumerate(columns)}

        samples = self.get_sample_list()

        mat = np.zeros((len(samples), len(columns)))
        mat[:] = np.nan

        for i, sample in enumerate(samples):
            for location, (unmethylated, methylated) in self.counts[sample].items():
                if dtype == 'methylated':
                    value = methylated
                elif dtype == 'unmethylated':
                    value = unmethylated
                elif dtype == 'beta':
                    value = methylated / (methylated + unmethylated)
                else:
                    raise ValueError
                mat[i, [column_to_index[location]]] = value

        return pd.DataFrame(mat, index=samples, columns=pd.MultiIndex.from_tuples(columns))

    def check_integrity(self):
        if len(self.sites) == 0 or len(self.counts) == 0:
            print(self)
            raise ValueError('The count matrix contains no data, verify if the input data was empty or filtered to stringently')


    def get_bulk_column(self, samples, location):

        total_un, total_met = 0, 0
        betas = []
        n_samples = 0
        for sample in samples:
            unmethylated, methylated = self.get_without_init((sample, location))
            total_un += unmethylated
            total_met += methylated

            if methylated + unmethylated > 0:
                n_samples += 1
                betas.append(methylated / (methylated + unmethylated))

        empty = (total_met+total_un) == 0
        return [ total_un, total_met, np.nan if empty else total_met/(total_un+total_met), np.var(betas) if len(betas) else np.nan, n_samples]


    def get_bulk_frame(self, dtype='pd', use_multi=True):
        """
        Get pandas dataframe containing the selected columns


        Returns:
            df(pd.DataFrame) : Dataframe containing the selected column, rows are locations,

        """
        self.check_integrity()
        # Fix columns
        columns = list(sorted(self.sites))
        # Create column to index mapping:
        column_to_index = {c: i for i, c in enumerate(columns)}

        samples = self.get_sample_list()

        mat = np.zeros((len(columns), 5))
        mat[:] = np.nan


        if use_multi and (self.threads is not None and self.threads>1):
            with Pool(self.threads) as workers:
                for index,column in enumerate(
                                                    workers.imap( get_bulk_vector,
                                                    ( (self, samples, location)
                                                    for index, location in enumerate(columns) ), chunksize=100_000)):
                    mat[index, :] =  column
        else:
            for index, location in enumerate(columns):
                mat[index, :] = self.get_bulk_column(samples, location)

        if dtype == 'pd':
            return pd.DataFrame(mat, index=pd.MultiIndex.from_tuples(columns),
                                columns=('unmethylated', 'methylated', 'beta', 'variance', 'n_samples'))
        elif dtype == 'np':
            return mat
        else:
            raise ValueError('dtype should be pd or np')


def methylation_dict_to_location_values(methylation_per_location_per_cell: dict, select_samples=None)->tuple:
    """
    Convert a dictionary
    { location -> cell -> [0,0] }
    into
    { contig : [ locations (list) ] }
    { contig : [ values (list) ] }
    """
    write_locations = defaultdict(list) # contig -> locations
    write_values = defaultdict(dict) # contig -> location -> value

    for location, cell_info_for_location in methylation_per_location_per_cell.items():
        # Calculate beta value:
        unmet = 0
        met = 0

        for cell, (c_unmet, c_met) in cell_info_for_location.items():
            if select_samples is not None and not cell in select_samples:
                continue
            unmet+=c_unmet
            met+=c_met

        support = unmet+met
        if support == 0:
            continue
        contig = location[0]
        position = location[1]

        write_locations[contig].append(position)
        write_values[contig][position] = met / support

    return write_locations, write_values



def twolist():
    return [0,0]

def defdict():
    return defaultdict(twolist)


def met_unmet_dict_to_betas(methylation_per_cell_per_cpg: dict, bin_size=None) -> dict:
    """
    Convert dictionary of count form to beta form:

    cell -> location -> [unmet, met]

    to

    cell -> location -> beta
    """
    export_table = defaultdict(dict) #location->cell->beta
    for (contig, start), data_per_cell in methylation_per_cell_per_cpg.items():
         for cell,(met,unmet) in data_per_cell.items():
                if type(start)==int and bin_size is not None:
                    export_table[cell][contig, start, start+bin_size] = met/ (unmet+met)
                else:
                    export_table[cell][contig, start] = met/ (unmet+met)
    return export_table


def extract_cpgs(bam,
                 contig,
                 fragment_class,
                 molecule_class,
                 start = None,
                 end = None,
                 fetch_start = None,
                 fetch_end = None,
                 context='Z',
                 stranded=False,
                 mirror_cpg = False,
                 allelic=False,
                 select_samples=None,
                 pool_alias = None,
                 reference_path = None,
                 methylation_consensus_kwargs= {},
                 bin_size=None):

    methylation_per_cell_per_cpg = defaultdict(defdict) # location -> cell -> [0,0]

    taps = TAPS()
    with pysam.AlignmentFile(bam) as al,\
         pysam.FastaFile((get_reference_path_from_bam(bam) if reference_path is None else reference_path)) as reference:

            for molecule in MoleculeIterator(
                al,
                fragment_class=fragment_class,
                molecule_class=molecule_class,
                molecule_class_args={
                     'reference':reference,
                     'taps':taps,
                     'taps_strand':'R',

                     'methylation_consensus_kwargs':methylation_consensus_kwargs,
                },
                fragment_class_args={},
                contig = contig,
                start=fetch_start,
                end=fetch_end
            ):
                if allelic:
                    allele =  molecule.allele

                if select_samples is not None and not molecule.sample in select_samples:
                    continue

                for (cnt, pos), call in molecule.methylation_call_dict.items():

                    if (start is not None and pos<start) or (end is not None and pos>=end):
                        continue

                    ctx = call['context']
                    if ctx.upper()!=context:
                        continue

                    if mirror_cpg and context=='Z' and not molecule.strand:
                        pos-=1

                    if pool_alias:
                        location_key = pool_alias
                    else:
                        if bin_size is not None:
                            location_key = [cnt, int(bin_size*int(pos/bin_size))]
                        else:
                            location_key = [cnt,pos]
                        if allelic:
                            location_key += [allele]

                        if stranded:
                            location_key += [molecule.strand]

                    methylation_per_cell_per_cpg[tuple(location_key)][molecule.sample][int(ctx.isupper())]+=1

    return methylation_per_cell_per_cpg
