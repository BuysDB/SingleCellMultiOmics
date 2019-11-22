import pandas as pd
import numpy as np
import glob
import math
import re
import sys
import multiprocessing


def downsampleRow(args):
    row, targetSum = args
    currentCount = row.sum()
    downsampledRow = row.copy()
    while currentCount > targetSum and currentCount != 0:
        possible = downsampledRow[(downsampledRow > 0)]
        desiredTossCount = int(currentCount - targetSum)
        probabilities = [p / currentCount for p in possible]
        for indexToLower in np.random.choice(
                possible.index, max(0, desiredTossCount),
                replace=True, p=probabilities):
            if downsampledRow[indexToLower] > 0:
                downsampledRow[indexToLower] -= 1

        currentCount = downsampledRow.sum()
    return downsampledRow

# downsample_to = sample to this amount of counts per column
# min_feature_abundance = remove all rows which have less than these counts


def downsampleDataFrame(df, downsample_to, min_feature_abundance=50):
    pool = multiprocessing.Pool(8)
    try:
        df = df.loc[:, df.sum() > downsample_to]
        df = df.loc[df.sum(1) > min_feature_abundance, :]
        subset = df.transpose()
        dfDownsampled = subset.copy()
        for idx, drow in enumerate(
            pool.map(
                downsampleRow, [
                (row, downsample_to) for i, row in subset.iterrows()])):
            dfDownsampled.iloc[idx, :] = drow
    except Exception as e:
        print(e)
    pool.close()
    return dfDownsampled.transpose()
