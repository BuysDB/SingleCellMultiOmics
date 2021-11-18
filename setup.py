#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))
version = {}
with open(os.path.join(_here, 'singlecellmultiomics', 'version.py')) as f:
    exec(f.read(), version)

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.md')) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

setup(

    name='singlecellmultiomics',
    version=version['__version__'],
    description='Tools to deal with one or more measurements from single cells',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Buys de Barbanson',
    author_email='b.barbanson@hubrecht.eu',
    url='https://github.com/BuysDB/SingleCellMultiOmics',
    download_url = 'https://github.com/BuysDB/SingleCellMultiOmics/archive/v0.1.9.tar.gz',

    license='MIT',
    packages=['singlecellmultiomics',

        'singlecellmultiomics.alleleTools',
        'singlecellmultiomics.bamProcessing',
        'singlecellmultiomics.barcodeFileParser',
        'singlecellmultiomics.countTableProcessing',
        'singlecellmultiomics.features',
        'singlecellmultiomics.fragment',
        'singlecellmultiomics.fastqProcessing',
        'singlecellmultiomics.fastaProcessing',
        'singlecellmultiomics.libraryDetection',
        'singlecellmultiomics.libraryProcessing',
        'singlecellmultiomics.modularDemultiplexer',
        'singlecellmultiomics.molecule',
        'singlecellmultiomics.methylation',
        'singlecellmultiomics.pyutils',
        'singlecellmultiomics.variants',
        'singlecellmultiomics.tags',
        'singlecellmultiomics.statistic',
        'singlecellmultiomics.tagtools',
        'singlecellmultiomics.universalBamTagger',
        'singlecellmultiomics.utils',
        'singlecellmultiomics.modularDemultiplexer.demultiplexModules'
        ],


    scripts=[
        # Demultiplexing
        'singlecellmultiomics/modularDemultiplexer/demux.py',

        # Fasta
        'singlecellmultiomics/fastaProcessing/fastaMaskVariants.py',
        'singlecellmultiomics/fastaProcessing/createMappabilityIndex.py',
        'singlecellmultiomics/fastaProcessing/fastaCreateDict.py',

        # Tagging
        'singlecellmultiomics/universalBamTagger/universalBamTagger.py',
        'singlecellmultiomics/universalBamTagger/bamtagmultiome_multi.py',
        'singlecellmultiomics/universalBamTagger/bamtagmultiome.py',
        'singlecellmultiomics/universalBamTagger/tapsTagger.py',
        'singlecellmultiomics/universalBamTagger/tapsTabulator.py',
        'singlecellmultiomics/universalBamTagger/4SUtagger.py',

        # Methylation:
        'singlecellmultiomics/methylation/methylationtab_to_bed.py',
        'singlecellmultiomics/methylation/bam_to_methylation_bw.py',

        # Bam processing:
        'singlecellmultiomics/bamProcessing/bamTabulator.py',
        'singlecellmultiomics/bamProcessing/bamToNucleosomePositions.py',
        'singlecellmultiomics/bamProcessing/bamAnalyzeCutDistances.py',
        'singlecellmultiomics/bamProcessing/bamSplitByTag.py',
        'singlecellmultiomics/bamProcessing/bamReadGroupFormat.py',
        'singlecellmultiomics/bamProcessing/bamSetReadGroup.py',
        'singlecellmultiomics/bamProcessing/bamToCountTable.py',
        'singlecellmultiomics/bamProcessing/bamCopyNumber.py',
        'singlecellmultiomics/bamProcessing/bamExtractSamples.py',
        'singlecellmultiomics/bamProcessing/scmoConvert.py',
        'singlecellmultiomics/bamProcessing/bamToMethylationAndCopyNumber.py',
        'singlecellmultiomics/bamProcessing/bamDinucleotideDistribution.py',
        'singlecellmultiomics/bamProcessing/bamToMethylationCalls.py',
        'singlecellmultiomics/bamProcessing/bamMappingRate.py',
        'singlecellmultiomics/bamProcessing/bamFilter.py',
        'singlecellmultiomics/bamProcessing/bamPlotRTstats.py',
        'singlecellmultiomics/bamProcessing/bamPlateVisualisation.py',
        'singlecellmultiomics/bamProcessing/bamFeatureDensityVisualisation.py',
        'singlecellmultiomics/bamProcessing/bamMethylationCutDistance.py',
        'singlecellmultiomics/bamProcessing/bamDuprate.py',
        'singlecellmultiomics/bamProcessing/bamExtractRandomPrimerStats.py',
        'singlecellmultiomics/bamProcessing/bamToRNACounts.py',
        'singlecellmultiomics/bamProcessing/bamToBigWig.py',
        'singlecellmultiomics/bamProcessing/structureTensor.py',
        'singlecellmultiomics/bamProcessing/bamCompare.py',
        'singlecellmultiomics/bamProcessing/variantStats.py',
        'singlecellmultiomics/bamProcessing/bamExtractVariants.py',
        'singlecellmultiomics/bamProcessing/bamMatchGATKBQSRReport.py',
        'singlecellmultiomics/bamProcessing/bamBinCounts.py',
        'singlecellmultiomics/bamProcessing/bamExtractNearMolecules.py',
        'singlecellmultiomics/bamProcessing/split_double_BAM.py',
        'singlecellmultiomics/bamProcessing/bamOverseq.py',
        'singlecellmultiomics/bamProcessing/plotRegion.py',
        'singlecellmultiomics/bamProcessing/soft_clip_filter.py',
        'singlecellmultiomics/utils/base_call_covariates.py',
        'singlecellmultiomics/bamProcessing/estimateTapsConversionEfficiency.py',

        # Library processing:
        'singlecellmultiomics/libraryProcessing/libraryStatistics.py',
        'singlecellmultiomics/libraryDetection/archivestats.py',
        'singlecellmultiomics/alleleTools/heterozygousSNPedit.py',

        # Feature conversion:
        'singlecellmultiomics/features/exonGTFtoIntronGTF.py',

        # Variants:
        'singlecellmultiomics/variants/postProcessVariants.py',
        'singlecellmultiomics/variants/vcfFilterAlleleFreq.py',
        'singlecellmultiomics/variants/vcfMutProfiler.py',
        'singlecellmultiomics/variants/plotCovariates.py',

        # Trimming:
        'singlecellmultiomics/fastqProcessing/trim_vasa.py',

        'singlecellmultiomics/utils/bigWigDiff.py',
        # Utility: (SGE wrapper)
        'singlecellmultiomics/utils/submission.py',
        'singlecellmultiomics/utils/ftp_upload.py',

        #Worfklow
        'singlecellmultiomics/snakemake_workflows/scmo_workflow.py',
        'singlecellmultiomics/snakemake_workflows/_general/sge_wrapper.py',
        'singlecellmultiomics/snakemake_workflows/_general/slurm_wrapper.py',

        # Facs:
        'singlecellmultiomics/FACS/trajectory.py'

        ],

  install_requires=[
       'pysam>=0.15.3','numpy>=1.16.5','pandas>=0.25.0','colorama','pyBigWig',
       'cutadapt>=2.9',
       'pysamiterators>=1.9','more-itertools','matplotlib','tabulate',
       'wheel','setuptools>=40.8.0','scikit-learn>=0.21.3','seaborn>=0.11.0', 'statsmodels', 'cached_property',
       'biopython>=1.71','pytest>=5.0.0','pytest-runner','snakemake>=5.8.1','lxml',
       'statsmodels' #,'tensorflow>=1.14.0'
   ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],

    include_package_data=True,

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
         'License :: OSI Approved :: MIT License'
        ]
)
