import pandas as pd

def dataframe_to_wig(df: pd.DataFrame, wig_path: str, span: int = 1, stepper: str = "variableStep",
                     graphType: str = "points", offset: int = 1):
    """
    Args:
        df: pandas.DataFrame to export to wig, expected format: multi_index (chrom,position)
            with columns: value,
        wig_path(str) : write WIG file to this location
        graphType: bar/points/line
        stepper:
        span:
    Returns:

    """
    with open(wig_path, 'w') as out:
        for chrom, chrom_data in df.sort_index().groupby(level=0):
            # Write the header:
            out.write(f'{stepper} chrom={chrom} graphType={graphType} span={span}\n')
            for k, row in chrom_data.iterrows():
                if len(k) == 2:
                    (chrom, pos) = k
                elif len(k) == 3:
                    (chrom, pos, _) = k
                else:
                    raise ValueError('Supply a dataframe with a multi index in the format (chrom, pos) ')
                out.write(f"{pos+offset}\t{row.iloc[0] if 'value' not in row else row['value']}\n")
