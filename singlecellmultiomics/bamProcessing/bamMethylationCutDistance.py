#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.rcParams['figure.dpi'] = 160
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import multiprocessing
from singlecellmultiomics.bamProcessing.bamBinCounts import generate_commands, count_methylation_binned
import argparse
from colorama import Fore, Style
from singlecellmultiomics.utils import dataframe_to_wig
from singlecellmultiomics.methylation import MethylationCountMatrix
from singlecellmultiomics.bamProcessing.bamFunctions import get_reference_from_pysam_alignmentFile
from colorama import Fore,Style
from collections import defaultdict, Counter
from multiprocessing import Pool
from datetime import datetime
import pysam
from singlecellmultiomics.bamProcessing import get_contig_sizes, get_contig_size
from singlecellmultiomics.bamProcessing.bamBinCounts import generate_commands, read_counts


def sample_dict():
    return defaultdict(Counter)



def methylation_to_cut_histogram(args):
    (alignments_path, bin_size, max_fragment_size, \
     contig, start, end, \
     min_mq, alt_spans, key_tags, dedup, kwargs) = args


    distance_methylation = defaultdict(sample_dict) # sample - > distance -> context(ZzHhXx) : obs
    max_dist = 1000

    # Define which reads we want to count:
    known =  set()
    if 'known' in kwargs and kwargs['known'] is not None:
        # Only ban the very specific TAPS conversions:
        try:
            with pysam.VariantFile(kwargs['known']) as variants:
                for record in variants.fetch(contig, start, end):
                    if record.ref=='C' and 'T' in record.alts:
                        known.add( record.pos)
                    if record.ref=='G' and 'A' in record.alts:
                        known.add(record.pos)
        except ValueError:
            # This happends on contigs not present in the vcf
            pass

    p = 0

    start_time = datetime.now()
    with pysam.AlignmentFile(alignments_path, threads=4) as alignments:
        # Obtain size of selected contig:
        contig_size = get_contig_size(alignments, contig)
        if contig_size is None:
            raise ValueError('Unknown contig')

        # Determine where we start looking for fragments:
        f_start = max(0, start - max_fragment_size)
        f_end = min(end + max_fragment_size, contig_size)

        for p, read in enumerate(alignments.fetch(contig=contig, start=f_start,
                                                  stop=f_end)):




            if p%50==0 and 'maxtime' in kwargs and kwargs['maxtime'] is not None:
                if (datetime.now() - start_time).total_seconds() > kwargs['maxtime']:
                    print(f'Gave up on {contig}:{start}-{end}')

                    break

            if not read_counts(read, min_mq=min_mq, dedup=dedup):
                continue


            tags = dict(read.tags)
            for i, (qpos, methylation_pos) in enumerate(read.get_aligned_pairs(matches_only=True)):

                # Don't count sites outside the selected bounds
                if methylation_pos < start or methylation_pos >= end:
                    continue

                call = tags['XM'][i]
                if call=='.':
                    continue

                sample = read.get_tag('SM')


                distance = abs(read.get_tag('DS') - methylation_pos)
                if distance>max_dist:
                    continue

                distance_methylation[sample][(read.is_read1, read.is_reverse, distance)][call] +=1

    return distance_methylation


threads = None

def get_distance_methylation(bam_path,
                                 bp_per_job: int,
                                 min_mapping_qual: int = None,
                                 skip_contigs: set = None,
                                 known_variants: str = None,
                                 maxtime: int = None,
                                 head: int=None,
                                 threads: int = None,
                                **kwargs
                                 ):


    all_kwargs = {'known': known_variants,
            'maxtime': maxtime,
            'threads':threads
            }
    all_kwargs.update(kwargs)
    commands = generate_commands(
        alignments_path=bam_path,
        key_tags=None,
        max_fragment_size=0,
        dedup=True,
        head=head,
        bin_size=bp_per_job,
        bins_per_job= 1, min_mq=min_mapping_qual,
        kwargs=all_kwargs,
        skip_contigs=skip_contigs
    )


    distance_methylation = defaultdict(sample_dict) # sample - > distance -> context(ZzHhXx) : obs

    with Pool(threads) as workers:

        for result in workers.imap_unordered(methylation_to_cut_histogram, commands):
            for sample, data_for_sample in result.items():
                for distance, context_obs in data_for_sample.items():
                    distance_methylation[sample][distance] += context_obs
    return distance_methylation


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract methylation levels relative to cut site (DS tag) from bam file""")

    argparser.add_argument('bamfile', metavar='bamfile', type=str)
    argparser.add_argument('-bp_per_job', default=5_000_000, type=int, help='Amount of basepairs to be processed per thread per chunk')
    argparser.add_argument('-threads', default=None, type=int, help='Amount of threads to use for counting, None to use the amount of available threads')

    fi = argparser.add_argument_group("Filters")
    fi.add_argument('-min_mapping_qual', default=40, type=int)
    fi.add_argument('-head', default=None, type=int,help='Process the first n bins')
    fi.add_argument('-skip_contigs', type=str, help='Comma separated contigs to skip', default='MT,chrM')
    fi.add_argument('-known_variants',
                           help='VCF file with known variants, will be not taken into account as methylated/unmethylated',
                           type=str)

    og = argparser.add_argument_group("Output")
    og.add_argument('-prefix', default='distance_calls', type=str, help='Prefix for output files')

    args = argparser.parse_args()



    print('Obtaining counts ', end="")
    r = get_distance_methylation(bam_path = args.bamfile,
                                 bp_per_job = args.bp_per_job,
                                 known_variants = args.known_variants,
                                 skip_contigs = args.skip_contigs.split(','),
                                 min_mapping_qual=args.min_mapping_qual,
                                 head = args.head,
                                 threads=args.threads,
    )
    print(f" [ {Fore.GREEN}OK{Style.RESET_ALL} ] ")


    for ctx in 'zhx':

        beta = {}
        met = {}
        un = {}
        for sample, sample_data in r.items():
            beta[sample] = {}
            met[sample] = {}
            un[sample] = {}
            for distance, contexts in sample_data.items():
                if ctx in contexts or ctx.upper() in contexts:
                    beta[sample][distance] = contexts[ctx.upper()]/(contexts[ctx.upper()]+contexts[ctx])
                    met[sample][distance] = contexts[ctx.upper()]
                    un[sample][distance] = contexts[ctx]

        pd.DataFrame(beta).sort_index().T.sort_index().to_csv(f'{args.prefix}_beta_{ctx}.csv')
        pd.DataFrame(beta).sort_index().T.sort_index().to_csv(f'{args.prefix}_beta_{ctx}.pickle.gz')
        pd.DataFrame(met).sort_index().T.sort_index().to_csv(f'{args.prefix}_counts_{ctx.upper()}.csv')
        pd.DataFrame(met).sort_index().T.sort_index().to_csv(f'{args.prefix}_counts_{ctx.upper()}.pickle.gz')
        pd.DataFrame(un).sort_index().T.sort_index().to_csv(f'{args.prefix}_counts_{ctx}.csv')
        pd.DataFrame(un).sort_index().T.sort_index().to_csv(f'{args.prefix}_counts_{ctx}.pickle.gz')

        # Make plots
        beta = {}
        met = {}
        un = {}
        for sample, sample_data in r.items():
            beta[sample] = {}
            met[sample] = {}
            un[sample] = {}
            for distance, contexts in sample_data.items():
                if distance[-1] > 500 or distance[-1] < 4: # Clip in sane region
                    continue
                if ctx in contexts or ctx.upper() in contexts:
                    beta[sample][distance] = contexts[ctx.upper()] / (contexts[ctx.upper()] + contexts[ctx])
                    met[sample][distance] = contexts[ctx.upper()]
                    un[sample][distance] = contexts[ctx]

        beta = pd.DataFrame(beta).sort_index().T.sort_index()
        met = pd.DataFrame(met).sort_index().T.sort_index()
        un = pd.DataFrame(un).sort_index().T.sort_index()

        for mate in [True, False]:
            for strand in [True, False]:
                fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

                un[mate, strand].sum().rename('Unmethylated').plot(ax=ax1)
                met[mate, strand].sum().rename('Methylated').plot(ax=ax1)

                ax1.set_xlabel('distance to cut')
                ax1.set_ylabel('# molecules')
                ax1.legend()

                (met[mate, strand].sum() / (un[mate, strand].sum() + met[mate, strand].sum())).rename('Beta').plot(
                    ax=ax2)
                # ax2.set_ylim(0,0)

                sns.despine()
                ax1.set_title(f'Mate {"R1" if mate else "R2"}, strand:{"reverse" if strand else "forward"}')
                ax2.set_ylabel('Beta')
                plt.savefig(f'{args.prefix}_{ctx}_{"R1" if mate else "R2"}_{"reverse" if strand else "forward"}.png')
                plt.close('all')
