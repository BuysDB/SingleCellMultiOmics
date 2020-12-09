import pandas as pd
from singlecellmultiomics.bamProcessing import get_contig_sizes
from collections import defaultdict
import pysam
import pyBigWig
import numpy as np

def dataframe_to_wig(df: pd.DataFrame, wig_path: str, span: int = 1, stepper: str = "variableStep",
                     graphType: str = "points", offset: int = 1, extra_header=''):
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
            out.write(f'{stepper} chrom={chrom} graphType={graphType} span={span} {extra_header}\n')
            for k, row in chrom_data.iterrows():
                if len(k) == 2:
                    (chrom, pos) = k
                elif len(k) == 3:
                    (chrom, pos, _) = k
                else:
                    raise ValueError('Supply a dataframe with a multi index in the format (chrom, pos) ')
                out.write(f"{pos+offset}\t{row.iloc[0] if 'value' not in row else row['value']}\n")


def write_bigwig(write_locations:dict, write_values:dict,  source_bam:str, write_path:str,bin_size: int=None):
    """
    Write location / values to a bigwig file

    Args:
        write_locations(dict) : Locations to write to, grouped by contig, {contig:[location (int), ..]

        write_values(dict) : Value to write, grouped by contig {contig:[value (float), ..],

        source_bam(str): path to source bam file to extract contig ordering from

        write_path(str): write bigwig file here

        bin_size(int) : bin_size, set to None to use a bin size of 1

    """
    with pysam.AlignmentFile(source_bam) as alignments, pyBigWig.open(write_path,'w') as out:

        cs = get_contig_sizes(alignments)
        # Write header
        out.addHeader(list(zip(alignments.references, alignments.lengths)))

        for contig in alignments.references:
            if not contig in write_locations:
                continue

            contig_list = [contig]*len(write_locations[contig])
            start = sorted(list(write_locations[contig]))
            if bin_size is None:
                end = np.array(start)+1
            else:
                end = np.array(start)+bin_size

            values_to_write = np.array( [write_values[contig][p] for p in start] , dtype=np.float32)

            out.addEntries(
                contig_list, #Contig
                start, #Start
                ends= end, #end
                values= values_to_write)
