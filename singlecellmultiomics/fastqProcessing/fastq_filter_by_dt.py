#!/usr/bin/env python3
import argparse
import gzip
from more_itertools import chunked

def select_datatype_fastq_pe(r1path, r2path, r1path_out, r2path_out, dt):
    look_for = f'dt:{dt}'
    with gzip.open(r1path,'rt') as r1h,\
        gzip.open(r2path,'rt') as r2h,\
        gzip.open(r1path_out,'wt') as r1ho,\
        gzip.open(r2path_out,'wt') as r2ho:

        for r1,r2 in zip(chunked(r1h,4), chunked(r2h,4)):
            if look_for in r1[0]:
                r1ho.write(''.join(r1))
                r2ho.write(''.join(r2))
                
def select_datatype_fastq_se(r1path, r1path_out, dt):
    look_for = f'dt:{dt}'
    with gzip.open(r1path,'rt') as r1h,\
        gzip.open(r1path_out,'wt') as r1ho:

        for r1 in chunked(r1h,4):
            if look_for in r1[0]:
                r1ho.write(''.join(r1))

                

def run():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Split fastq file based on dt: tag"
    )

    argparser.add_argument(
        '-r1_in',
        type=str,required=True
        )
    argparser.add_argument(
        '-r2_in',
        type=str,required=False
        )
    argparser.add_argument(
        '-r1_out',help='R1 output',
        type=str,required=True
        )
    argparser.add_argument(
        '-r2_out',help='R2 output',
        type=str,required=False
        )
    
    argparser.add_argument(
        '-dt',
        help="Datatype, for example RNA or DamID",
        type=str,required=True
        )
    
    
    args = argparser.parse_args()
    if args.r2_in is not None:
        select_datatype_fastq_pe(args.r1_in, 
                          args.r2_in, 
                          args.r1_out, 
                          args.r2_out, 
                          args.dt)
    else:
         select_datatype_fastq_se(args.r1_in, 
                          args.r1_out, 
                          args.dt)
    print('Done')

if __name__ == '__main__':
    run()