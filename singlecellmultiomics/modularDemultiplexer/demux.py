#!/usr/bin/env python3
# -*- coding: utf-8 -*-
f"!!! PLEASE USE PYTHON 3.6 OR HIGHER !!!"
import os
import sys
import re
import gzip
from colorama import Fore
from colorama import Back
from colorama import Style
import argparse
from singlecellmultiomics.fastqProcessing.fastqHandle import FastqHandle
import singlecellmultiomics.barcodeFileParser.barcodeFileParser as barcodeFileParser
import singlecellmultiomics.libraryDetection.sequencingLibraryListing as sequencingLibraryListing
from singlecellmultiomics.modularDemultiplexer.demultiplexingStrategyLoader import DemultiplexingStrategyLoader
import glob
from colorama import init
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import NonMultiplexable

if __name__=='__main__':

	# Location of this script
	demuxer_location = os.path.dirname(os.path.realpath(__file__))

	init()
	import logging
	logging.getLogger().setLevel(logging.WARNING)
	#logging.getLogger().setLevel(logging.INFO)

	argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""%sMulti-library Single cell-demultiplexing%s, written by Buys de Barbanson, Hubrecht institute 2016-2018
										%sExample usage:%s

										List how the demultiplexer will interpret the data: (No demultiplexing is executed)
										%sdemultiplex.py ./*.fastq.gz%s
										Demultiplex a set of fastq files into separate files per cell:
										%sdemultiplex.py ./*.fastq.gz --y%s
										Demultiplex on the cluster:
										%sdemultiplex.py ./*.fastq.gz -submit demux.sh; sh demux.sh %s
		"""% (Style.BRIGHT, Style.RESET_ALL, Style.BRIGHT, Style.RESET_ALL, Style.DIM, Style.RESET_ALL, Style.DIM, Style.RESET_ALL, Style.DIM, Style.RESET_ALL))

	argparser.add_argument('fastqfiles', metavar='fastqfile', type=str, nargs='*', help='Input fastq files. For windows compatibility, if any of the file names contains a star the expression is evaluated as GLOB. Read names should be in illumina format or LIBRARY_R1.fastq.gz, LIBRARY_R2.fastq.gz, if one file is supplied and the name does not end on .fastq. or fastq.gz the file is interpreted as a list of files')
	argparser.add_argument('--y', help="Perform the demultiplexing (Accept current parameters)", action='store_true' )
	argparser.add_argument('-experimentName', '-e', help="Name of the experiment, set to uf to use the folder name", type=str, default='uf' )

	inputArgs = argparser.add_argument_group('Input', '')
	inputArgs.add_argument('-n', help="Only get the first n properly demultiplexable reads from a library", type=int)
	inputArgs.add_argument('-r', help="Only get the first n reads from a library, then demultiplex", type=int)
	inputArgs.add_argument('-slib', help="Assume all files belong to the same library, this flag supplies the name", type=str )
	inputArgs.add_argument('-replace', action='append', default=None, help="""Replace part of the library name by another string, [SEARCH,REPLACEMENT]. For example if you want to remove FOO from the library name use "-replace FOO," if you want to replace FOO by BAR use "-replace FOO,BAR" """)

	inputArgs.add_argument('-merge', default=None, help="""Merge libraries through multiple runs by selecting a merging method.
	Options:None, delimiter[position], examples, given the samplename 'library_a_b': -merge _ yields 'library', -merge _2 yields 'library_a'as identifier.""")
	inputArgs.add_argument('--ignore', action='store_true', help="Ignore non-demultiplexable read files")
	#inputArgs.add_argument('-bfsp', help="Barcode file searchpaths", type=str, default='/media/sf_data/references,/hpc/hub_oudenaarden/bdebarbanson/ref')
	outputArgs = argparser.add_argument_group('Output', '')
	argparser.add_argument('--norejects', help="Do not store rejected reads",  action='store_true')

	outputArgs.add_argument('-o', help="Output (cell) file directory, when not supplied the current directory/raw_demultiplexed is used", type=str, default='./raw_demultiplexed')
	outputArgs.add_argument('--scsepf', help="Every cell gets a separate FQ file", action='store_true' )


	bcArgs = argparser.add_argument_group('Barcode', '')
	bcArgs.add_argument('-hd', help="Hamming distance barcode expansion; accept cells with barcodes N distances away from the provided barcodes. Collisions are dealt with automatically. ", type=int, default=0)
	bcArgs.add_argument('--lbi', help="List barcodes being used for cell demultiplexing", action='store_true')

	bcArgs.add_argument('-barcodeDir', default=os.path.join(demuxer_location,'barcodes'), help="Directory from which to obtain the barcodes")



	indexArgs = argparser.add_argument_group('Sequencing indices', '')
	indexArgs.add_argument('--li', help="List sequencing indices.", action='store_true')
	indexArgs.add_argument('-indexDir', default=os.path.join(demuxer_location,'indices'), help="Directory from which to obtain the sequencing indices")
	indexArgs.add_argument('-si', help="Select only these sequencing indices -si CGATGT,TTAGGC")
	indexArgs.add_argument('-hdi', help="Hamming distance INDEX sequence expansion, the hamming distance used for resolving the sequencing INDEX. For cell barcode hamming distance see -hd", type=int, default=1)

	fragArgs = argparser.add_argument_group('Fragment configuration', '')
	#fragArgs.add_argument('--rc1', help="First mate is reverse complemented", action='store_true')
	#fragArgs.add_argument('--rc2', help="Second mate is reverse complemented", action='store_true')
	fragArgs.add_argument('--se', help="Allow single end reads",  action='store_true')

	techArgs = argparser.add_argument_group('Technical', '')
	#techArgs.add_argument('-t', help="Amount of demultiplexing threads used" , type=int, default=8)
	techArgs.add_argument('-fh', help="When demultiplexing to mutliple cell files in multiple threads, the amount of opened files can exceed the limit imposed by your operating system. The amount of open handles per thread is kept below this parameter to prevent this from happening.", default=500, type=int)
	techArgs.add_argument('-dsize', help="Amount of reads used to determine barcode type" , type=int, default=10000)

	argparser.add_argument('-use',default=None, help='use these demultplexing strategies, comma separate to select multiple. For example for cellseq 2 data with 6 basepair umi: -use CS2C8U6 , for combined mspji and Celseq2: MSPJIC8U3,CS2C8U6 if nothing is specified, the best scoring method is selected' )


	argparser.add_argument('-ignoreMethods', help='Do not try to load these methods', default="scarsMiSeq"  )

	argparser.add_argument('-maxAutoDetectMethods','-mxa', help='When --use is not specified, how many methods can the demultiplexer choose at the same time? This loosely corresponds to the amount of measurements you made in a single cell', default=1, type=int  )
	argparser.add_argument('-minAutoDetectPct','-mia', help='When --use is not specified, what is the lowest percentage yield required to select a demultplexing strategy', default=2, type=float  )
	args = argparser.parse_args()
	verbosity = 1

	ignoreMethods = args.ignoreMethods.split(',')

	if len(set(args.fastqfiles))!=len(args.fastqfiles):
		print(f'{Fore.RED}{Style.BRIGHT}Some fastq files are supplied multiple times! Pruning those!{Style.RESET_ALL}')
		args.fastqfiles = set(args.fastqfiles)

	if len(args.fastqfiles)==1:
		if not args.fastqfiles[0].endswith('.gz') and not args.fastqfiles[0].endswith('.fastq')  and not args.fastqfiles[0].endswith('.fq'):
			# File list:
			print('Input is interpreted as a list of files..')
			with open(args.fastqfiles[0]) as f:
				fqFiles = []
				for line in f:
					fqFiles.append( line.strip() )
			args.fastqfiles = fqFiles

	# Load barcodes
	barcodeParser = barcodeFileParser.BarcodeParser(hammingDistanceExpansion=args.hd, barcodeDirectory=args.barcodeDir)




	## Setup the index parser
	indexFileAlias=None # let the multiplex methods decide which index file to use
	if args.si: ## the user defines the sequencing indices
		useSequencingIndices = args.si.split(',')
		print(f"{Style.BRIGHT} Only these sequencing indices will be kept: {Style.RESET_ALL}")
		for sequencingIndex in useSequencingIndices:
			print(f"{Fore.GREEN}{sequencingIndex}{Style.RESET_ALL}")

		indexParser = barcodeFileParser.BarcodeParser()
		indexFileAlias='user'
		for index,sequencingIndex in enumerate(useSequencingIndices):
			indexParser.addBarcode(
				index = str(index),
				barcodeFileAlias='user',
				barcode=sequencingIndex,
				hammingDistance=0,
				originBarcode=None)
		# Perform expansion:
		indexParser.expand(args.hdi, alias='user')

	else: # the sequencing indices are automatically detected
		indexParser =  barcodeFileParser.BarcodeParser(hammingDistanceExpansion=args.hdi, barcodeDirectory=args.indexDir)

	if args.lbi:
		barcodeParser.list(showBarcodes=None)
	if args.li:
		indexParser.list(showBarcodes=None)

	#Load the demultiplexing strategies
	dmx = DemultiplexingStrategyLoader(barcodeParser=barcodeParser, indexParser=indexParser,ignoreMethods=ignoreMethods, indexFileAlias=indexFileAlias)
	dmx.list()

	if len(args.fastqfiles)==0:
		print(f'{Fore.RED}No files supplied, exitting.{Style.RESET_ALL}')
		exit()

	print(f"\n{Style.BRIGHT}Detected libraries:{Style.RESET_ALL}")
	libraries = sequencingLibraryListing.SequencingLibraryLister().detect(args.fastqfiles, args)

	# Detect the libraries:
	if args.use is None:
		if len(libraries)==0:
			raise ValueError('No libraries found')

		print(f"\n{Style.BRIGHT}Demultiplexing method Autodetect results{Style.RESET_ALL}")
		# Run autodetect
		processedReadPairs, strategyYieldsForAllLibraries = dmx.detectLibYields(libraries, testReads=args.dsize, maxAutoDetectMethods=args.maxAutoDetectMethods, minAutoDetectPct=args.minAutoDetectPct)

	print(f"\n{Style.BRIGHT}Demultiplexing:{Style.RESET_ALL}")
	for library in libraries:
		if args.use is None:
			processedReadPairs = strategyYieldsForAllLibraries[library]['processedReadPairs']
			strategyYieldForLibrary =  strategyYieldsForAllLibraries[library]['strategyYields']
			selectedStrategies = dmx.selectedStrategiesBasedOnYield(processedReadPairs, strategyYieldForLibrary, maxAutoDetectMethods = args.maxAutoDetectMethods, minAutoDetectPct=args.minAutoDetectPct)
			selectedStrategies = dmx.getSelectedStrategiesFromStringList(selectedStrategies)
		else:
			selectedStrategies = dmx.getSelectedStrategiesFromStringList(args.use.split(','))

		print(f'Library {library} will be demultiplexed using:')
		for stra in selectedStrategies:
			print(f'\t{Fore.GREEN}{str(stra)}{Style.RESET_ALL}')
		if len(selectedStrategies)==0:
			print(f'{Fore.RED}NONE! The library will not be demultiplexed!{Style.RESET_ALL}')

		if not args.y:
			#with open(args.submit, 'w') as f:

			filesForLib = []
			for lane in libraries[library]:
				for R1R2 in libraries[library][lane]:
					for p in libraries[library][lane][R1R2]:
						filesForLib.append( p )
			arguments = " ".join([x for x in sys.argv if x!='--dry' and not '--y' in x and not '-submit' in x and not '.fastq' in x and not '.fq' in x]) + " --y"

			print(f"\n{Style.BRIGHT}--y not supplied, execute the command below to run demultiplexing on the cluster:{Style.RESET_ALL}")
			print( '/hpc/hub_oudenaarden/bdebarbanson/internalTools/submission.py' + f' -y --nenv -time 50 -t 1 -m 8 -N NDMX%s "source /hpc/hub_oudenaarden/bdebarbanson/virtualEnvironments/py36/bin/activate; %s -use {",".join([x.shortName for x in selectedStrategies])}"\n' % (library, '%s %s'  % ( arguments, " ".join(filesForLib)) ))

		if args.y:
			targetDir = f'{args.o}/{library}'
			if not os.path.exists(targetDir):
				os.makedirs(targetDir)
			handle = FastqHandle(f'{args.o}/{library}/demultiplexed' , True, single_cell=args.scsepf, maxHandles=args.fh)

			rejectHandle = FastqHandle(f'{args.o}/{library}/rejects' , True)

			processedReadPairsForThisLib = 0

			for lane, readPairs in libraries[library].items():
				if args.n and processedReadPairsForThisLib>=args.n:
					break
				for readPair in readPairs:
					pass
				for readPairIdx,_ in enumerate(readPairs[readPair]):
					files = [ readPairs[readPair][readPairIdx] for readPair in readPairs ]
					processedReadPairs,strategyYields = dmx.demultiplex( files , strategies=selectedStrategies, targetFile=handle, rejectHandle=rejectHandle,
					library=library, maxReadPairs=None if args.n is None else (args.n-processedReadPairsForThisLib))
					processedReadPairsForThisLib += processedReadPairs
					if args.n and processedReadPairsForThisLib>=args.n:
						break
			handle.close()
			rejectHandle.close()
