from setuptools import setup
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.rst')) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, 'README.rst'), encoding='utf-8') as f:
        long_description = f.read()

version = {}
with open(os.path.join(_here, '\singlecelmultiomics', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name='singlecelmultiomics',
    version=version['__version__'],
    description=('Tools to deal with one or more measurements from single cells'),
    long_description=long_description,
    author='Buys de Barbanson',
    author_email='b.barbanson@hubrecht.eu',
    url='https://github.com/bast/somepackage',
    license='MPL-2.0',
    packages=['singlecelmultiomics'],

  install_requires=[
       'pysam==0.15.0'
   ],
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.6'],
)
