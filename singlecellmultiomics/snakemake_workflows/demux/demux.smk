from singlecellmultiomics.libraryDetection.sequencingLibraryListing import SequencingLibraryLister
from glob import glob
import collections

# This code detects which libraries are present in the current folder:
l = SequencingLibraryLister()
LIBRARIES = l.detect(glob('*.fastq.gz'), merge='_')
# Flatten:
fastq_per_lib = collections.defaultdict(list)
for lib,lane_dict in LIBRARIES.items():
    for lane,read_dict in lane_dict.items():
        fastq_per_lib[lib] += read_dict['R1']
        fastq_per_lib[lib] += read_dict['R2']
libraries =  list( fastq_per_lib.keys() )
###

def get_fastq_file_list(wildcards):
    # Obtain a list of fastq files associated to wildcards.library
    global libraries
    return sorted( fastq_per_lib[wildcards.library] )

def get_target_demux_list():
    global libraries
    targets = []
    for lib in libraries:
        targets.append('processed/' + lib + "/demultiplexedR1.fastq.gz" )
        targets.append('processed/' + lib + "/demultiplexedR2.fastq.gz" )
    return targets


rule all:
    input:
        get_target_demux_list()

rule demux:
    input:
        fastqfiles = get_fastq_file_list
    output:
        "processed/{library}/demultiplexedR1.fastq.gz",
        "processed/{library}/demultiplexedR2.fastq.gz"
    shell:
        "demux.py -merge _ {input.fastqfiles} -o processed --y -n 10000"
