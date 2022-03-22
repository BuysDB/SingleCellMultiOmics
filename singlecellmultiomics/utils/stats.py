import statsmodels.api as sm
import numpy as np
from scipy import stats
from copy import copy
import pandas as pd
import statsmodels.formula.api as smf
from singlecellmultiomics.utils import pool_wrapper
from multiprocessing import Pool
from scipy.ndimage import gaussian_filter
import pyBigWig
from typing import Callable, Union

def calculate_nested_f_statistic(small_model, big_model):
    # From https://stackoverflow.com/a/60769343/2858160
    """Given two fitted GLMs, the larger of which contains the parameter space of the smaller, return the F Stat and P value corresponding to the larger model adding explanatory power"""
    addtl_params = big_model.df_model - small_model.df_model
    f_stat = (small_model.deviance - big_model.deviance) / (addtl_params * big_model.scale)
    df_numerator = addtl_params
    # use fitted values to obtain n_obs from model object:
    df_denom = (big_model.fittedvalues.shape[0] - big_model.df_model)
    p_value = stats.f.sf(f_stat, df_numerator, df_denom)
    return (f_stat, p_value)


def gaussian_2d(df:pd.DataFrame, sigmas:tuple, **kwargs) -> pd.DataFrame:
    """
    Smooth a pd.DataFrame using a gaussian filter (scipy.ndimage),
    kwargs are passed to the gaussian_filter function


    """
    df = pd.DataFrame(
        gaussian_filter(df, sigmas, **kwargs),
        index=df.index,
        columns=df.columns)
    return df


def vectorize_bw(path: str, bin_size: int, selected_contigs:list ,sample=None) -> pd.Series:
    """Vectorize (numpy) a bigwig file

    Args:
        path: path to bigwig file
        bin_size : bin size
        selected_contigs: list of contigs to take into account
        sample: name of output series, tuple for MultiIndex
     """

    bin_labels=[]
    densities= []
    with pyBigWig.open(path) as handle:
        for contig, size in handle.chroms().items():
            if not contig in selected_contigs:
                continue
            nbins = int(size/bin_size)+1
            density = np.zeros( nbins)
            for bin_index, bin_start in enumerate(range(0,size,bin_size)):
                bin_end = bin_start + bin_size
                bin_end=  min(bin_end,size-1)
                try:
                    density[bin_index] = np.array( handle.intervals(contig,bin_start, bin_end) )[:,2].mean()
                except IndexError:
                    pass

                bin_labels.append((contig,bin_start,bin_end))
            densities.append(density)
    # Convert to series:
    return  pd.Series(np.concatenate(densities),index=bin_labels,name=sample)

def vectorize_bw_multi(paths: list, selected_contigs:list, formatsample: Callable[ [str,], Union[tuple,str,int] ], bin_size:int, n_threads=None):

    commands = [
        (
            vectorize_bw,
            {
                'sample' : formatsample(path),
                'path' : path,
                'bin_size': bin_size,
                'selected_contigs':selected_contigs
            }
        )
        for path in paths
    ]
    with Pool(n_threads) as workers:
        densities=[]
        for ser in workers.imap(
            pool_wrapper,
            commands
        ):
            densities.append(ser)

    return pd.concat(densities,axis=1).T


def _GLM_cluster_de_test_single_gene(gene, cuts_frame, clusters):

    """
    Calculate if a gene varies due to batch effects or is significantly changing between clusters
    """

    data = copy(cuts_frame[[gene]])+1
    data.columns = ['ncuts']
    data['plate'] = [x.split('_')[0] for x in data.index]
    data['cluster'] = clusters
    data['n_total_cuts'] = cuts_frame.sum(1)


    fam = sm.families.Poisson()
    try:
        model = smf.glm("ncuts ~ 1 + plate + cluster", data= data,
                 # cov_struct=ind,
                offset=np.log(data['n_total_cuts']),
                      family=fam).fit()

        null_model = smf.glm(f"ncuts ~ 1 + plate", data= data,
                 # cov_struct=ind,
                offset=np.log(data['n_total_cuts']),
                      family=fam).fit()
    except Exception as e:
        if 'estimation infeasible.'  in str(e) or 'PerfectSeparationError' in str(e) :
            return None
        else:
            raise

    coeff = pd.DataFrame( {'model_std_err':model.bse,
               'model_coefficients':model.params,
               'null_std_err':null_model.bse,
               'null_coefficients':null_model.params,

              })

    return [gene, *calculate_nested_f_statistic(null_model,model), coeff, model, null_model]


def GLM_cluster_de_test(cuts_frame, clusters):

    """
    Calculate if a gene varies due to batch effects or is significantly changing between clusters
    """

    table = []
    for gene in cuts_frame.columns:

        r = _GLM_cluster_de_test_single_gene(gene, 1+cuts_frame, clusters)

        if r is None:
            continue

        (gene, f_score, p_value, coeff, model, null_model) = r
        table.append([gene,f_score,p_value, coeff['model_std_err']['cluster'], coeff['model_coefficients']['cluster']])

    return pd.DataFrame(table, columns=['gene','f_score','p_value','cluster_stderr','cluster_coeff'])


def GLM_cluster_de_test_multi(df, y, n_processes=None):

    """
    Calculate if a gene varies due to batch effects or is significantly changing between clusters, multiprocessed

    Args:
        df(pd.DataFrame) : dataframe of cuts,  rows are cells, columns are loci, values are cut-counts (not normalised)

        y(np.array) : target vector / clusters

        n_processes(int) : Amount of processes to use
    """
    table = []
    with Pool(n_processes) as workers:

        for r in workers.imap_unordered(pool_wrapper,((

            _GLM_cluster_de_test_single_gene
        ,{
            'gene':gene,
            'cuts_frame':df,
            'clusters':y
        }
        )
        for gene in df
        ), chunksize=200):
            if r is None:
                continue

            (gene, f_score, p_value, coeff, model, null_model) = r
            table.append([gene,f_score,p_value, coeff['model_std_err']['cluster'], coeff['model_coefficients']['cluster']])


    return pd.DataFrame(table, columns=['gene','f_score','p_value','cluster_stderr','cluster_coeff'])
