#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np


class MethylationCountMatrix:

    def __init__(self, counts: dict = None):
        # Sample->(contig,bin_start,bin_end)-> [methylated_counts, unmethylated]
        self.counts = {} if counts is None else counts
        self.sites = set()

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
        sample, location = key
        if sample not in self.counts or location not in self.counts[sample]:
            return [0, 0]
        return self.counts[sample][location]

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

        for location, row in self.get_bulk_frame().iterrows():
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

    def get_bulk_frame(self, dtype='pd'):
        """
        Get pandas dataframe containing the selected columns


        Returns:
            df(pd.DataFrame) : Dataframe containing the selected column, rows are locations,

        """
        # Fix columns
        columns = list(sorted(self.sites))
        # Create column to index mapping:
        column_to_index = {c: i for i, c in enumerate(columns)}

        samples = self.get_sample_list()

        mat = np.zeros((len(columns), 5))
        mat[:] = np.nan

        for index, location in enumerate(columns):

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

            mat[index, 0:2] = (total_un, total_met)
            if mat[index, 0:2].sum() > 0:
                mat[index, 2] = total_met / (total_un + total_met)

            if len(betas) > 0:
                mat[index, 3] = np.var(betas)

            mat[index, 4] = n_samples
        if dtype == 'pd':
            return pd.DataFrame(mat, index=pd.MultiIndex.from_tuples(columns),
                                columns=('unmethylated', 'methylated', 'beta', 'variance', 'n_samples'))
        elif dtype == 'np':
            return mat
        else:
            raise ValueError('dtype should be pd or np')






