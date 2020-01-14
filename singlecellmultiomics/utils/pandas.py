import pandas as pd
import numpy as np
import seaborn as sns

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
