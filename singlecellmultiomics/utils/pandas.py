import pandas as pd
import numpy as np
import seaborn as sns
import warnings

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
