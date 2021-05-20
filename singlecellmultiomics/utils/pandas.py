import pandas as pd
import numpy as np
import seaborn as sns
import warnings
from scipy.stats import zscore


def z_score(df, axis):
    """Calculate z-score of pandas dataframe"""
    return pd.DataFrame(zscore(df,axis), index=df.index, columns=df.columns)

def linear_scale(df,vmin,vmax,fmin=None,fmax=None):
    """
    Scale the values in the dataframe using One-dimensional linear interpolation.
    """
    if fmin is None:
        fmin = df.min().min()

    if fmax is None:
        fmax = df.max().max()

    return pd.DataFrame(
        np.interp(df,
            (fmin,fmax),
            (vmin, vmax)),
        index=df.index,
        columns=df.columns
        )

def rolling_smooth(df: pd.DataFrame, window:int, center=True, axis=0):
    """Rolling window mean smoothing of pandas dataframe, clips of the missing values"""
    if center:
        offset = int(np.ceil(window/2))
        if axis==0:
            return df.rolling(window,center=center,axis=axis).mean().iloc[offset:-offset]
        else:
            return df.rolling(window,center=center,axis=axis).mean().iloc[:,offset:-offset]
    else:
        if axis==0:
            return df.rolling(window,center=center,axis=axis).mean().iloc[window:]
        else:
            return df.rolling(window,center=center,axis=axis).mean().iloc[:,window:]

def createRowColorDataFrame( discreteStatesDataFrame, nanColor =(0,0,0), predeterminedColorMapping={} ):

    """ Create color dataframe for use with seaborn clustermap

    Args:
        discreteStatesDataFrame (pd.DataFrame) : Dataframe containing the data to convert to colors, like:  pd.DataFrame( [['A','x'],['A','y']],index=['A','B'], columns=['First', 'Second'] )

        nanColor(tuple) : Color for records having an NAN

        predeterminedColorMapping(dict) : Supply class colors here (optional)

    Returns:
        discreteColorMatrix (pd.DataFrame) : Dataframe to pass to seaborn clustermap row_colors, or col_colors

        luts (dict) : class->color mapping
    """
    # Should look like:
    # discreteStatesDataFrame = pd.DataFrame( [['A','x'],['A','y']],index=['A','B'], columns=['First', 'Second'] )
    colorMatrix = []
    luts = {}
    for column in discreteStatesDataFrame:
        states = [x for x in discreteStatesDataFrame[column].unique() if not pd.isnull(x)]
        undeterminedColorStates = [x for x in discreteStatesDataFrame[column].unique() if not pd.isnull(x) and not x in predeterminedColorMapping]

        cols = sns.color_palette('hls',len(undeterminedColorStates))
        #lut = { i:sns.color_palette('bright').jet(x) for i,x in zip(states, np.linspace(0,1,len(states)) )}
        lut = { state:cols[i] for i,state in enumerate(undeterminedColorStates) }
        lut.update({key:value for key,value in predeterminedColorMapping.items() if key in states})
        lut[np.nan] = nanColor
        colorMatrix.append( [ nanColor if pd.isnull(x) else lut[x] for x in  discreteStatesDataFrame[column] ] )
        luts[column] = lut
    discreteColorMatrix = pd.DataFrame(colorMatrix, index=discreteStatesDataFrame.columns, columns=discreteStatesDataFrame.index ).transpose()
    return discreteColorMatrix, luts

def interpolate(series, target):
    """
    Interpolate the given pd.series at the coordinates given by target (np.array)
    the index of the series is the x coordinate, the values xp
    """
    # prune missing data:
    series = series[~pd.isnull(series)]

    return  pd.Series(
        #interpolate:
        np.interp(target, series.index,  series.values),
        # re use existing series name:
        name=series.name,
        # use target as new index:
        index=target
    )


def tordist(x1: float, x2: float, wrap_dist: float ) -> float:
    """Calculate the toroidial distance between two scalars

    Args:
        x1(float) : first datapoint
        x2(float) : second datapoint
        wrap_dist(float) : wrapping distance (highest value), values higher than this will wrap around to zero

    Returns:
        distance(float) : toroidial distance between x1 and x2, wrapping around wrap_dist
    """
    dx = abs(x2 - x1)
    if dx>wrap_dist*0.5:
        return wrap_dist-dx
    else:
        return dx

def tor_resample(x: np.array, y: np.array, window_radius: float, max_tp: float,n:int=100) -> pd.Series:
    """ Toroidal resample a set of coordinates x,y, where x is a set of timepoints into a new set of coordinates from zero to max_tp with n steps. Uses a sliding mean."""
    interp = {}
    s = pd.Series(y,index=x)

    warnings.simplefilter("ignore")
    for tp in np.linspace(0,max_tp, n):

        selected_points = np.array([( tordist(x,tp,max_tp) <= window_radius) for x,y in s.items()])

        q = s[selected_points]
        mean = np.nanmean(q)
        interp[tp] = mean
        interp[tp-max_tp] = mean
        interp[tp+max_tp] = mean

    resampled = pd.Series(interp).sort_index()
    return resampled.loc[0:max_tp]
