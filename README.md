Single cell multi omics is a set of tools to deal with multiple measurements from the same cel.
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
