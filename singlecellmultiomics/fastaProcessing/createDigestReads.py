#!/usr/bin/env python3
import singlecellmultiomics.utils
import pysam
import collections
import numpy as np
import re
import more_itertools
import gzip
import argparse
import uuid
import os

if self.__name__ == '__main__':
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""In silico digest genome""")
    argparser.add_argument('fasta', metavar='fastafile', type=str, help='Fasta file to digest')

    argparser.add_argument('-minlen', type=int, default=20, help='minimum read length')
    argparser.add_argument('-maxlen', type=int, default=69, help='maximum read length')
    argparser.add_argument('-head', type=int)
    argparser.add_argument('-contigs', type=str, default=None, help='Comma separated contigs to run analysis for, when None specified all contigs are used')
    args = argparser.parse_args()

    r1_read_length = args.maxlen
    minlen = args.minlen

    ref = pysam.FastaFile(args.fasta)

    def seq_to_fastq(seq, header):
        return f'@{header}\n{seq}\n+\n{"E"*len(seq)}\n'

    selected_contigs = None if args.contigs is None else set(args.contigs.split(','))

    fastq_path = f'simulated_nlaIII_single_{r1_read_length}_{uuid.uuid4()}.fastq.gz'

    outbam = f'simulated_nlaIII_single_{r1_read_length}.bam'




    with gzip.open(fastq_path,'wt') as out:
        processed = 0


        for contig, contig_len in zip(ref.references, ref.lengths):
            if selected_contigs is not None and not contig in selected_contigs:
                continue
            print(f'{contig}\t{contig_len}')
            seq = ref.fetch(contig).upper()
            frag_locations = np.diff( [m.start() for m in re.finditer('CATG', seq)])
            # Generate fragments:
            for i,(start, end_excluding_catg) in enumerate(more_itertools.windowed(  (m.start() for m in re.finditer('CATG', seq)), 2 )):
                if end_excluding_catg is None or args.head is not None and i>=(args.head-1):
                    continue
                end = end_excluding_catg+4


                forward_read = seq[start:end][:r1_read_length]
                if len(forward_read)>=minlen:
                    out.write(
                        seq_to_fastq(forward_read, f'DR:fwd;ct:{contig};ds:{start};qn:{len(forward_read)}')
                    )
                reverse_read = singlecellmultiomics.utils.reverse_complement(seq[start:end])[:r1_read_length]
                if len(reverse_read)>=minlen:
                    out.write(
                        seq_to_fastq(reverse_read, f'DR:rev;ct:{contig};ds:{end_excluding_catg};qn:{len(reverse_read)}')
                    )

            processed += contig_len

    # Map the fastq file:
    os.system(f"bwa mem -t 4 {args.fasta} {fastq_path} | samtools view -b - > ./{outbam}.unsorted.bam; samtools sort -T ./temp_sort -@ 4 ./{outbam}.unsorted.bam > ./{outbam}.unfinished.bam & mv ./{outbam}.unfinished.bam ./{outbam} & samtools index {outbam} & rm {fastq_path}")

    # Create site dictionary:
    sites = {}  #site-> wrongly_assinged_to, correctly_assigned, lost

    ref_reads = pysam.AlignmentFile(outbam)
    # Perform site detection:
    for read in ref_reads:

        if read.is_supplementary:
            continue
        # Check if there is reads mapping to this site

        kv = {kv.split(':')[0]:kv.split(':')[1] for kv in read.query_name.split(';')}

        site_pos = int(kv['ds'])
        site_chrom = kv['ct']
        strand = (kv['DR']=='rev')

        read_site =   read.reference_end-4 if read.is_reverse else read.reference_start
        read_contig = read.reference_name

        key = (site_chrom, site_pos, strand)
        key_mapped = (read_contig, read_site, read.is_reverse)

        if not key in sites:
            sites[key] = {'correct':0, 'wrong_gain':0, 'lost':0}

        if read_site == site_pos and site_chrom==read_contig and read.is_reverse==strand: # Correct assignment
            sites[key]['correct'] +=1
        else:
            # Assign a missing count to the origin site:
            sites[key]['lost'] +=1

            # Assing a wrong annotation to the target site:
            if not key_mapped in sites:
                sites[key_mapped] =  {'correct':0, 'wrong_gain':0, 'lost':0}
            sites[key_mapped]['wrong_gain'] += 1


    outtab = f'simulated_nlaIII_single_{r1_read_length}.mappability.stats.tsv'
    outtabsafe = f'simulated_nlaIII_single_{r1_read_length}.mappability.safe.tsv'

    with open(outtab,'w') as o, open(outtabsafe,'w') as of:
        for (contig,pos,strand), measured in sites.items():
            o.write(f'{contig}\t{pos}\t{"+-"[strand]}\t{measured["correct"]}\t{measured["lost"]}\t{measured["wrong_gain"]}\n')

            if measured['wrong_gain']==0 and measured['lost']==0 and measured['correct']==1:
                of.write(f'{contig}\t{pos}\t{"+-"[strand]}\n')
