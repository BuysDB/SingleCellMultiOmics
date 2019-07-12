[![Build Status](https://travis-ci.com/BuysDB/SingleCellMultiOmics.svg?branch=master)](https://travis-ci.com/BuysDB/SingleCellMultiOmics)

## Single cell multi omics
Single cell multi omics is a set of tools to deal with multiple measurements from the same cell.

# Installation
```
git clone https://github.com/BuysDB/HandleLimiter
git clone https://github.com/BuysDB/pysamiterators
git clone https://github.com/BuysDB/SingleCellMultiOmics
pip install -e ./HandleLimiter  -e ./pysamiterators -e ./SingleCellMultiOmics
```

# Usage

The following tools are available:

demux.py:
This tool demultiplexes raw fastq files, adapters, molecule and cell information are removed from the fastq records and encoded into the read name including their base-calling qualities.
Additional stored inforation includes:
- Assigned and raw illumina index
- Library
- Demultiplexing strategy used for demultiplexing (What kind of data the read is derived from)
- Assigned barcode index

The next step is usually to trim and then map the demultiplexed reads.

 For RNA seq data aligned to a transcriptome the step after this is to run featureCounts.

The mapped reads are encoded in a BAM file. This BAM file still contains the encoded data and this has to be decoded in order to get a useful BAM file.
universalBamTagger.py
1) Recodes the original read names and extracts all information previously encoded by the demultiplexer.
2) Adds allele information. (A VCF file is required for this)
3) Supports multiple protocols: RNA:CELSEQ1 and CELSEQ2 (with 8 and 6bp UMI), methylation digest sequencing:SC MSPJI ,  lineage tracing:SCARTRACE, DNA digest sequencing: NLAIII, histone modification sequencing: scCHIC.
4) Assigns reads to molecules to allow for deduplication, adds duplication BAM flag
5) Assigns read groups



# Examples:

Demultiplex all fastq.gz files in the current directory using NLAIII barcodes
```
demux.py *.fastq.gz -use NLAIII384C8U3 --y
````

Demultiplex only the specified sequencing index (GTTTGA), and everything 1 hamming distance away from GTTTGA  :
```
demux.py -si GTTTGA *.gz --y --hdi 1
```

### scCHIC
For every fragment in input.bam find scCHIC seq fragments and deduplicate these. Fragments with the same cell barcode, umi, library and strand and starting within a range of 5 bp from each other are assigned as duplicate. The mnase cut site location is expected to be between the first base (Usually an A) this A is part of the sequencing adapter, and the second base (Usually a T). The cut site location is recorded into the DS tag. When alleles are specified using -alleles, the molecule assignment is split up by allele, this means that if two fragments map to the same location and share the same umi, but contain SNPs which indicate differing alleles, the reads are not assigned to the same molecule. For every fragment the ligation sequence is recorded into the RZ tag.
```
universalBamTagger.py --chic --ftag -moleculeRadius 5  -o tagged.bam input.bam
 ```

[Complete scCHIC data processing instructions from FastQ to count table here](https://github.com/BuysDB/SingleCellMultiOmics/wiki/scCHIC-data-processing)

### NlaIII
For every fragment in input.bam find NLAIII seq fragments and deduplicate these. Fragments with the same cell barcode, umi, library and strand are assigned as duplicate. The NlaIII cut site location is recorded into the DS tag. When alleles are specified using -alleles, the molecule assignment is split up by allele, this means that if two fragments map to the same location and share the same UMI, but contain SNPs which indicate differing alleles, the reads are not assigned to the same molecule. For every fragment the sequenced part of the NlaIII cut site sequence is recorded into the RZ tag, this is usually CATG, but is allowed to be shifted 1 base to ATG. In the NlaIII protocol a reverse transcription (RT) is used, generally capturing more reverse transcription reactions will yield a more accurate molecule consensus sequence. For every fragment which support the molecule the reverse transcription reaction is recorded by storing the location of the random primer used for RT and the sequence of the random primer.
```
universalBamTagger.py --nla --ftag  -o tagged.bam input.bam
 ```
 [Complete NlaIII data processing instructions from FastQ to count table here](https://github.com/BuysDB/SingleCellMultiOmics/wiki/NLA-III-data-processing)

### Molecule iteration in Python
Assigned molecules can be iterated over to use in your Python code.
Every molecule is a list of fragments which are assigned to the molecule.

```
with pysam.AlignmentFile('example.bam') as f:
    for molecule in MoleculeIterator(f):
        print(molecule)
```
[More detailed examples here](https://github.com/BuysDB/SingleCellMultiOmics/wiki/Molecule-iteration).

### Plate visualisation

Show relative abundance of reads and unique molecules across 384 well plate.
```
bamPlateVisualisation.py tagged.bam -o ./plate_plots
```
Creates the folder ./plate_plots containing  
```
raw_reads_[LIBRARY_TYPE]_[LIBRARY_NAME].png # Distribution of total reads
usable_reads_[LIBRARY_TYPE]_[LIBRARY_NAME].png # Distribution of reads which can be assigned to a molecule
unique_reads_[LIBRARY_TYPE]_[LIBRARY_NAME].png # Distribution of unique reads
```


Create a contig by sample matrix and divide counts when reads are multimapping. (Used for counting transcriptome mapped reads)
```
bamToCountTable.py -featureTags chrom -sampleTags SM --divideMultimapping --dedup ./tagged/STAR_mappedAligned.sortedByCoord.out.bam -o transcriptome_counts.csv
```

Obtain sample, chromosome, restrictionsite, read start, and read end from test.bam file:
```
bamFileTabulator.py -featureTags SM,reference_name,DS,reference_start,reference_end test.bam
```
List all available tags:
```
bamFileTabulator.py test.bam
```

You can additionaly use any of these pysam read attributes:
```
aend
alen
aligned_pairs
bin
blocks
cigarstring
cigartuples
compare
flag
header
inferred_length
is_duplicate
is_paired
is_proper_pair
is_qcfail
is_read1
is_read2
is_reverse
is_secondary
is_supplementary
is_unmapped
isize
mapping_quality
mapq
mate_is_reverse
mate_is_unmapped
mpos
mrnm
next_reference_id
next_reference_name
next_reference_start
opt
overlap
pnext
pos
positions
qend
qlen
qname
qqual
qstart
qual
query
query_alignment_end
query_alignment_length
query_alignment_qualities
query_alignment_sequence
query_alignment_start
query_length
query_name
query_qualities
query_sequence
reference_end
reference_id
reference_length
reference_name
reference_start
rlen
rname
rnext
seq
template_length
tid
tlen
```
