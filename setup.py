from setuptools import setup
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.md')) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

version = {}
with open(os.path.join(_here, 'singlecellmultiomics', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name='singlecelmultiomics',
    version=version['__version__'],
    description=('Tools to deal with one or more measurements from single cells'),
    long_description=long_description,
    author='Buys de Barbanson',
    author_email='b.barbanson@hubrecht.eu',
    url='https://github.com/BuysDB/SingleCellMultiOmics',
    license='MPL-2.0',
    packages=['singlecellmultiomics'],
    scripts=[
        'singlecellmultiomics/modularDemultiplexer/demux.py',

        'singlecellmultiomics/universalBamTagger/universalBamTagger.py',

        # Bam processing:
        'singlecellmultiomics/bamProcessing/bamTabulator.py',
        'singlecellmultiomics/bamProcessing/bamToCountTable.py',
        'singlecellmultiomics/bamProcessing/bamFilter.py',
        'singlecellmultiomics/bamProcessing/bamPlateVisualisation.py',

        # Library processing:
        'singlecellmultiomics/libraryProcessing/libraryStatistics.py'

        ],
  install_requires=[
       'pysam'
   ],
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.6'],
)
