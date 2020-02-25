#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pysam
import argparse
import sys
import collections
import pysam
import pandas as pd
from uuid import uuid4


def tables_gen(path,verbose=False):
    with open(path) as i:
        iterable=iter(i)
        header = next(i)
        table_data=None
        n_columns = 0
        n_rows = 0
        formatting=None
        table_name=None
        table_header=None
        for line in iterable:
            line = line.strip()

            if line.startswith('#'):
                if table_data is not None:
                    yield table_generator(formatting, table_name, table_header, table_data, n_columns,n_rows)
                formatting = line.split(':')
                n_columns,n_rows = int(formatting[2]), int(formatting[3])
                line = next(i).strip()
                table_name = line.split(':')
                table_header =  next(i).strip().split()
                if verbose:
                    print('FMT',formatting)
                    print('table_name',table_name)
                    print('table_header',table_header)
                table_data = []
                continue
            # fill the table:
            parts = line.strip().split()
            if len(parts)==n_columns:
                table_data.append(parts)
        if formatting is not None and table_name is not None and table_header is not None:
            yield table_generator(formatting, table_name, table_header, table_data, n_columns,n_rows)

def table_generator(formatting, table_name, table_header, table_data, n_columns,n_rows):
    """
    Generates pandas dataframes from GATK formatted tables
    """
    df = pd.DataFrame(table_data)
    try:
        df.columns = table_header[-n_columns:]
    except Exception as e:
        #print(e)
        pass
    return [df,
            {
                'formatting':formatting,
                'table_name':table_name,
                'table_header':table_header,
                'n_columns':n_columns,
                'n_rows':n_rows

            }]

def match_bam_with_GATK_BQSR_report(input_bam_path, output_bam_path, bqsr_report_path ):

    # Identify all read groups present in the bam:
    with pysam.AlignmentFile(input_bam_path) as a:
        read_groups_in_bam = set(rg['ID'] for rg in a.header['RG'] )

    print( f'{len(read_groups_in_bam)} read groups present in bam file')

    print('Now reading BQSR report...')
    # Find which read groups are missing in the bqsr report file:
    total_missing= set()
    for i,(table, meta) in enumerate(tables_gen(path=bqsr_report_path)):
        if 'ReadGroup' in table.columns:
            read_groups_in_gatk_table = set( table['ReadGroup'].unique() )
            missing = read_groups_in_bam.difference(read_groups_in_gatk_table)
            total_missing.update(missing)
        del table
        del meta

    rg_file_path = f'./{str( uuid4() )}.rg.txt'
    total_present = read_groups_in_bam.difference(total_missing)
    print(f'{len(total_missing)} read groups are missing in the report. Keeping {len(total_present)}')
    # Write bam file with only selected read groups
    with open(rg_file_path,'w') as out:
        for rg in total_present:
            out.write(f'{rg}\n')

    print('Extracting alignments ...')

    os.system( f'samtools view -h {input_bam_path} -R {rg_file_path}  -O BAM -o {output_bam_path} -@ 8' )
    os.remove(rg_file_path)
    print('All done.')


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract samples/cells from bam file, writes to multiple bam files if multiple groups are supplied
    """)
    argparser.add_argument('bamfile', metavar='bam_input_file', type=str)
    argparser.add_argument('bqsr', metavar='bqsr_file', type=str)
    argparser.add_argument('-o', type=str, help='output path bam', required=True)
    args = argparser.parse_args()

    match_bam_with_GATK_BQSR_report(args.bamfile,args.o, args.bqsr )
