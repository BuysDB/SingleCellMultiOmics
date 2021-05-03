import statsmodels.api as sm
import numpy as np
from scipy import stats
from copy import copy
import pandas as pd
import statsmodels.formula.api as smf

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


def _GLM_cluster_de_test_single_gene(gene, cuts_frame, clusters):

    """
    Calculate if a gene varies due to batch effects or is significantly changing between clusters
    """

    data = copy(cuts_frame[[gene]])
    data.columns = ['ncuts']
    data['plate'] = [x.split('_')[0] for x in data.index]
    data['cluster'] = clusters
    data['n_total_cuts'] = cuts_frame.sum(1)

    fam = sm.families.Poisson()

    model = smf.glm("ncuts ~ 1 + plate + cluster", data= data,
             # cov_struct=ind,
            offset=np.log(data['n_total_cuts']),
                  family=fam).fit()

    null_model = smf.glm(f"ncuts ~ 1 + plate", data= data,
             # cov_struct=ind,
            offset=np.log(data['n_total_cuts']),
                  family=fam).fit()

    return [gene, *calculate_nested_f_statistic(null_model,model)]

def GLM_cluster_de_test(cuts_frame, clusters):

    """
    Calculate if a gene varies due to batch effects or is significantly changing between clusters
    """
    return pd.DataFrame(
        [
        _GLM_cluster_de_test_single_gene(gene, cuts_frame, clusters)
        for gene in cuts_frame.columns
        ],columns=['gene','fstat','pval'])
