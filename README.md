## Single cell multi omics
Single cell multi omics is a set of tools to deal with multiple measurements from the same cell.

# Installation
```
pip install https://github.com/BuysDB/HandleLimiter/archive/master.zip
pip install https://github.com/BuysDB/pysamiterators/archive/master.zip
pip install https://github.com/BuysDB/SingleCellMultiOmics/archive/master.zip
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
3) Supports multiple protocols: RNA:CELSEQ1 and CELSEQ2 (with 8 and 6bp UMI), methylation digest sequencing:SC MSPJI ,  lineage tracing:SCARTRACE, DNA digest sequencing: NLAIII, histone modification sequencing: CHIC.
4) Assigns reads to molecules to allow for deduplication, adds duplication BAM flag
5) Assigns read groups



# Examples:

Demultiplex all fastq.gz files in the current directory using NLAIII barcodes
```
demux.py *.fastq.gz -use NLAIII384C8U3 --y
````

For every fragment in input.bam find if it is a valid CHIC seq fragment and deduplicate. Fragments with the same barcode and umi and starting 5 bp from each other are assigned as duplicate.
```universalBamTagger.py --chic --ftag -moleculeRadius 5  -o tagged.bam input.bam ```

Create a contig by sample matrix and divide counts when reads are multimapping. (Used for counting transcriptome mapped reads)
```bamToCountTable.py -featureTags chrom -sampleTags SM --divideMultimapping --dedup ./tagged/STAR_mappedAligned.sortedByCoord.out.bam -o transcriptome_counts.csv```

Obtain SM and DS tags from tagged bam file ( Sample and restriction site)
```bamFileTabulator.py SM DS test.bam```


Create a table
```bamFileToCountTable.py ```

 	Added and renamed bam processing tools 	3 hours ago
bamFilter.py
