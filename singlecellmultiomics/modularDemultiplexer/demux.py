#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pkg_resources
import logging
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import NonMultiplexable
from colorama import init
from singlecellmultiomics.modularDemultiplexer.demultiplexingStrategyLoader import DemultiplexingStrategyLoader
import singlecellmultiomics.libraryDetection.sequencingLibraryListing as sequencingLibraryListing
import singlecellmultiomics.barcodeFileParser.barcodeFileParser as barcodeFileParser
from singlecellmultiomics.fastqProcessing.fastqHandle import FastqHandle
import argparse
from colorama import Style
from colorama import Fore
import sys
import os
f"!!! PLEASE USE PYTHON 3.6 OR HIGHER !!!"

if __name__ == '__main__':

    # Location of this script
    demuxer_location = os.path.dirname(os.path.realpath(__file__))

    init()

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""%sMulti-library Single cell-demultiplexing%s, written by Buys de Barbanson, Hubrecht institute 2016-2018
										%sExample usage:%s

										List how the demultiplexer will interpret the data: (No demultiplexing is executed)
										%sdemultiplex.py ./*.fastq.gz%s
										Demultiplex a set of fastq files into separate files per cell:
										%sdemultiplex.py ./*.fastq.gz --y%s
										Demultiplex on the cluster:
										%sdemultiplex.py ./*.fastq.gz -submit demux.sh; sh demux.sh %s
		""" %
        (Style.BRIGHT,
         Style.RESET_ALL,
         Style.BRIGHT,
         Style.RESET_ALL,
         Style.DIM,
         Style.RESET_ALL,
         Style.DIM,
         Style.RESET_ALL,
         Style.DIM,
         Style.RESET_ALL))

    argparser.add_argument(
        'fastqfiles',
        metavar='fastqfile',
        type=str,
        nargs='*',
        help='Input fastq files. Read names should be in illumina format or LIBRARY_R1.fastq.gz, LIBRARY_R2.fastq.gz, if one file is supplied and the name does not end on .fastq. or fastq.gz the file is interpreted as a list of files')
    argparser.add_argument(
        '--y',
        help="Perform the demultiplexing (Accept current parameters)",
        action='store_true')
    argparser.add_argument(
        '-experimentName',
        '-e',
        help="Name of the experiment, set to uf to use the folder name",
        type=str,
        default='uf')

    inputArgs = argparser.add_argument_group('Input', '')
    inputArgs.add_argument(
        '-n',
        help="Only get the first n properly demultiplexable reads from a library",
        type=int)
    inputArgs.add_argument(
        '-r',
        help="Only get the first n reads from a library, then demultiplex",
        type=int)
    inputArgs.add_argument(
        '-slib',
        help="Assume all files belong to the same library, this flag supplies the name",
        type=str)
    inputArgs.add_argument(
        '-replace',
        action='append',
        default=None,
        help="""Replace part of the library name by another string, [SEARCH,REPLACEMENT]. For example if you want to remove FOO from the library name use "-replace FOO," if you want to replace FOO by BAR use "-replace FOO,BAR" """)

    inputArgs.add_argument(
        '-merge',
        default="_",
        help="""Merge libraries through multiple runs by selecting a merging method.
	Options:None, delimiter[position], examples, given the samplename 'library_a_b': -merge _ yields 'library', -merge _2 yields 'library_a'as identifier.""")
    inputArgs.add_argument(
        '--ignore',
        action='store_true',
        help="Ignore non-demultiplexable read files")
    #inputArgs.add_argument('-bfsp', help="Barcode file searchpaths", type=str, default='/media/sf_data/references,/hpc/hub_oudenaarden/bdebarbanson/ref')
    outputArgs = argparser.add_argument_group('Output', '')
    argparser.add_argument(
        '--norejects',
        help="Do not store rejected reads",
        action='store_true')

    outputArgs.add_argument(
        '-o',
        help="Output (cell) file directory, when not supplied the current directory/raw_demultiplexed is used",
        type=str,
        default='./raw_demultiplexed')
    outputArgs.add_argument(
        '--scsepf',
        help="Every cell gets a separate FQ file",
        action='store_true')
    outputArgs.add_argument(
        '-nextcmd',
        help="Execute this command when the demultiplexing is finished. When cluster submission is used this command is executed after all jobs are finished",
        type=str,
        default=None)

    bcArgs = argparser.add_argument_group('Barcode', '')
    bcArgs.add_argument(
        '-hd',
        help="Hamming distance barcode expansion; accept cells with barcodes N distances away from the provided barcodes. Collisions are dealt with automatically. ",
        type=int,
        default=0)
    bcArgs.add_argument(
        '--lbi',
        help="List barcodes being used for cell demultiplexing",
        action='store_true')
    bcArgs.add_argument(
        '-barcodeDir',
        default=pkg_resources.resource_filename(
            'singlecellmultiomics',
            'modularDemultiplexer/barcodes/'),
        help="Directory from which to obtain the barcodes, when nothing is supplied the package resources are used")

    indexArgs = argparser.add_argument_group('Sequencing indices', '')
    indexArgs.add_argument(
        '--li',
        help="List sequencing indices.",
        action='store_true')
    indexArgs.add_argument(
        '-indexDir',
        default=pkg_resources.resource_filename(
            'singlecellmultiomics',
            'modularDemultiplexer/indices/'),
        help="Directory from which to obtain the sequencing indices,  when nothing is supplied the package resources are used")
    indexArgs.add_argument(
        '-si', help="Select only these sequencing indices -si CGATGT,TTAGGC")
    indexArgs.add_argument(
        '-ifa',
        help="Index file alias, select from the list supplied when running --li",
        default='illumina_merged_ThruPlex48S_RP',
        type=str)
    indexArgs.add_argument(
        '-hdi',
        help="Hamming distance INDEX sequence expansion, the hamming distance used for resolving the sequencing INDEX. For cell barcode hamming distance see -hd",
        type=int,
        default=1)

    fragArgs = argparser.add_argument_group('Fragment configuration', '')
    #fragArgs.add_argument('--rc1', help="First mate is reverse complemented", action='store_true')
    #fragArgs.add_argument('--rc2', help="Second mate is reverse complemented", action='store_true')
    fragArgs.add_argument(
        '--se',
        help="Allow single end reads",
        action='store_true')

    techArgs = argparser.add_argument_group('Technical', '')
    techArgs.add_argument(
        '-g',
        help="group_id, don't use this yourself",
        type=int,
        default=None)
    techArgs.add_argument(
        '-fh',
        help="When demultiplexing to mutliple cell files in multiple threads, the amount of opened files can exceed the limit imposed by your operating system. The amount of open handles per thread is kept below this parameter to prevent this from happening.",
        default=500,
        type=int)
    techArgs.add_argument(
        '-dsize',
        help="Amount of reads used to determine barcode type",
        type=int,
        default=2000)
    techArgs.add_argument(
        '--nochunk',
        help="Do not run lanes in separate jobs",
        action='store_true')

    argparser.add_argument(
        '-use',
        default=None,
        help='use these demultplexing strategies, comma separate to select multiple. For example for cellseq 2 data with 6 basepair umi: -use CS2C8U6 , for combined mspji and Celseq2: MSPJIC8U3,CS2C8U6 if nothing is specified, the best scoring method is selected')

    argparser.add_argument(
        '-ignoreMethods',
        help='Do not try to load these methods',
        default="scarsMiSeq")

    argparser.add_argument(
        '-maxAutoDetectMethods',
        '-mxa',
        help='When --use is not specified, how many methods can the demultiplexer choose at the same time? This loosely corresponds to the amount of measurements you made in a single cell',
        default=1,
        type=int)
    argparser.add_argument(
        '-minAutoDetectPct',
        '-mia',
        help='When --use is not specified, what is the lowest percentage yield required to select a demultplexing strategy',
        default=2,
        type=float)
    args = argparser.parse_args()
    verbosity = 1

    ignoreMethods = args.ignoreMethods.split(',')
    if any(('*' in fq_file for fq_file in args.fastqfiles)):
        raise ValueError(
            "One or more of the fastq file paths contain a '*', these files cannot be interpreted, review your input files")

    if len(set(args.fastqfiles)) != len(args.fastqfiles):
        print(f'{Fore.RED}{Style.BRIGHT}Some fastq files are supplied multiple times! Pruning those!{Style.RESET_ALL}')
        args.fastqfiles = set(args.fastqfiles)

    if len(args.fastqfiles) == 1:
        if not args.fastqfiles[0].endswith('.gz') and not args.fastqfiles[0].endswith(
                '.fastq') and not args.fastqfiles[0].endswith('.fq'):
            # File list:
            print('Input is interpreted as a list of files..')
            with open(args.fastqfiles[0]) as f:
                fqFiles = []
                for line in f:
                    fqFiles.append(line.strip())
            args.fastqfiles = fqFiles

    # Sort the fastq files..
    args.fastqfiles = sorted(args.fastqfiles)

    # Load barcodes
    barcodeParser = barcodeFileParser.BarcodeParser(
        hammingDistanceExpansion=args.hd,
        barcodeDirectory=args.barcodeDir)

    # Setup the index parser
    indexFileAlias = args.ifa  # let the multiplex methods decide which index file to use

    if args.si:  # the user defines the sequencing indices
        useSequencingIndices = args.si.split(',')
        print(
            f"{Style.BRIGHT} Only these sequencing indices will be kept: {Style.RESET_ALL}")
        for sequencingIndex in useSequencingIndices:
            print(f"{Fore.GREEN}{sequencingIndex}{Style.RESET_ALL}")

        indexParser = barcodeFileParser.BarcodeParser()
        indexFileAlias = 'user'
        for index, sequencingIndex in enumerate(useSequencingIndices):
            indexParser.addBarcode(
                index=str(index),
                barcodeFileAlias='user',
                barcode=sequencingIndex,
                hammingDistance=0,
                originBarcode=None)
        # Perform expansion:
        indexParser.expand(args.hdi, alias='user')

    else:  # the sequencing indices are automatically detected
        indexParser = barcodeFileParser.BarcodeParser(
            hammingDistanceExpansion=args.hdi, barcodeDirectory=args.indexDir)

    if args.lbi:
        barcodeParser.list(showBarcodes=None)
    if args.li:
        indexParser.list(showBarcodes=None)

    # Load the demultiplexing strategies
    dmx = DemultiplexingStrategyLoader(
        barcodeParser=barcodeParser,
        indexParser=indexParser,
        ignoreMethods=ignoreMethods,
        indexFileAlias=indexFileAlias)
    dmx.list()

    if len(args.fastqfiles) == 0:
        print(f'{Fore.RED}No files supplied, exitting.{Style.RESET_ALL}')
        exit()

    print(f"\n{Style.BRIGHT}Detected libraries:{Style.RESET_ALL}")
    libraries = sequencingLibraryListing.SequencingLibraryLister().detect(
        args.fastqfiles, args=args)

    # Detect the libraries:
    if args.use is None:
        if len(libraries) == 0:
            raise ValueError('No libraries found')

        print(
            f"\n{Style.BRIGHT}Demultiplexing method Autodetect results{Style.RESET_ALL}")
        # Run autodetect
        processedReadPairs, strategyYieldsForAllLibraries = dmx.detectLibYields(
            libraries, testReads=args.dsize, maxAutoDetectMethods=args.maxAutoDetectMethods, minAutoDetectPct=args.minAutoDetectPct, verbose=True)

    print(f"\n{Style.BRIGHT}Demultiplexing:{Style.RESET_ALL}")
    for library in libraries:
        if args.use is None:
            processedReadPairs = strategyYieldsForAllLibraries[library]['processedReadPairs']
            strategyYieldForLibrary = strategyYieldsForAllLibraries[library]['strategyYields']
            selectedStrategies = dmx.selectedStrategiesBasedOnYield(
                processedReadPairs,
                strategyYieldForLibrary,
                maxAutoDetectMethods=args.maxAutoDetectMethods,
                minAutoDetectPct=args.minAutoDetectPct)
            selectedStrategies = dmx.getSelectedStrategiesFromStringList(
                selectedStrategies)
        else:
            selectedStrategies = dmx.getSelectedStrategiesFromStringList(
                args.use.split(','))

        print(f'Library {library} will be demultiplexed using:')
        for stra in selectedStrategies:
            print(f'\t{Fore.GREEN}{str(stra)}{Style.RESET_ALL}')
        if len(selectedStrategies) == 0:
            print(
                f'{Fore.RED}NONE! The library will not be demultiplexed! The used barcodes could not be detected automatically. Please supply the desired method using the -use flag or increase the -dsize parameter to use more reads for detecting the library type.{Style.RESET_ALL}')

        if not args.y:
            print(f"\n{Style.BRIGHT}--y not supplied, execute the command below to run demultiplexing on the cluster:{Style.RESET_ALL}")
            arguments = " ".join(
                [x for x in sys.argv if x != '--dry' and '--y' not in x and '-submit' not in x and '.fastq' not in x and '.fq' not in x]) + " --y"

            submit_in_chunks = (not args.scsepf and not args.nochunk)
            submitted_jobs = []
            filesForLib = []
            group_id = 0
            for lane in libraries[library]:
                files_to_submit = []
                for R1R2 in libraries[library][lane]:

                    for p in libraries[library][lane][R1R2]:
                        filesForLib.append(p)
                        files_to_submit.append(p)

                if submit_in_chunks:
                    job_name = f'DMX_{library}_{group_id}'
                    submitted_jobs.append(job_name)

                    print(
                        'submission.py' +
                        f' -y -sched auto --py36 -time 50 -t 1 -m 8 -N {job_name} "%s  -g {group_id} -use {",".join([x.shortName for x in selectedStrategies])}"\n' %
                        ('%s %s' %
                         (arguments,
                          " ".join(files_to_submit))))
                group_id += 1

            final_jobs = []
            if not submit_in_chunks:
                job_name = f'DMX_{library}'
                print(
                    'submission.py' +
                    f' -y --py36 -time 50 -t 1 -m 8 -sched auto -N {job_name} "%s -use {",".join([x.shortName for x in selectedStrategies])}"\n' %
                    ('%s %s' %
                     (arguments,
                      " ".join(filesForLib))))
                final_jobs.append(job_name)
            else:
                # we need a job which glues everything back together
                # f'{args.o}/{library}/{prefix}demultiplexed
                # f'{args.o}/{library}/{prefix}rejects
                cmds = [
                    f'cat {args.o}/{library}/*_TEMP_demultiplexedR1.fastq.gz  > {args.o}/{library}/demultiplexedR1.fastq.gz && rm {args.o}/{library}/*_TEMP_demultiplexedR1.fastq.gz',
                    f'cat {args.o}/{library}/*_TEMP_demultiplexedR2.fastq.gz  > {args.o}/{library}/demultiplexedR2.fastq.gz && rm {args.o}/{library}/*_TEMP_demultiplexedR2.fastq.gz',
                    f'cat {args.o}/{library}/*_TEMP_demultiplexing.log  > {args.o}/{library}/demultiplexing.log && rm {args.o}/{library}/*_TEMP_demultiplexing.log']
                if not args.norejects:
                    cmds += [
                        f'cat {args.o}/{library}/*_TEMP_rejectsR1.fastq.gz  > {args.o}/{library}/rejectsR1.fastq.gz && rm {args.o}/{library}/*_TEMP_rejectsR1.fastq.gz',
                        f'cat {args.o}/{library}/*_TEMP_rejectsR2.fastq.gz  > {args.o}/{library}/rejectsR2.fastq.gz && rm {args.o}/{library}/*_TEMP_rejectsR2.fastq.gz'
                    ]

                for i, cmd in enumerate(cmds):
                    job_name = f'glue_{library}_{i}'
                    print(
                        'submission.py' +
                        f' -y -sched auto --silent --py36 -time 4 -t 1 -m 2 -N "glue_{library}" "{cmd}" -hold {",".join(submitted_jobs)}')
                    final_jobs.append(job_name)

            # Execute last command if applicable
            if args.nextcmd is not None:
                print(
                    'submission.py' +
                    f' -y --silent -sched auto -time 1 -t 1 -m 1 -N "NEXT_{library}" "{cmd}" -hold {",".join(final_jobs)}')

        if args.y:

            targetDir = f'{args.o}/{library}'
            if not os.path.exists(targetDir):
                try:
                    os.makedirs(targetDir)
                except FileExistsError:
                    continue

            prefix = '' if args.g is None else f'{args.g}_TEMP_'
            handle = FastqHandle(
                f'{args.o}/{library}/{prefix}demultiplexed',
                True,
                single_cell=args.scsepf,
                maxHandles=args.fh)
            if args.norejects:
                rejectHandle = None
            else:
                rejectHandle = FastqHandle(
                    f'{args.o}/{library}/{prefix}rejects', True)
            """Set up statistic file"""

            log_location = os.path.abspath(
                f'{args.o}/{library}/{prefix}demultiplexing.log')
            log_handle = open(log_location, 'w')
            log_handle.write(" ".join(sys.argv) + '\n')
            log_handle.write(
                f'Demultiplexing operation started, writing to {args.o}/{library}\n')

            processedReadPairsForThisLib = 0
            for lane, readPairs in libraries[library].items():
                if args.n and processedReadPairsForThisLib >= args.n:
                    break
                for readPair in readPairs:
                    pass

                for readPairIdx, _ in enumerate(readPairs[readPair]):
                    files = [readPairs[readPair][readPairIdx]
                             for readPair in readPairs]
                    log_handle.write(
                        f"processing input files:\t{','.join(files)}\n")
                    processedReadPairs, strategyYields = dmx.demultiplex(files,
                                                                         strategies=selectedStrategies,
                                                                         targetFile=handle,
                                                                         rejectHandle=rejectHandle,
                                                                         log_handle=log_handle,
                                                                         library=library,
                                                                         maxReadPairs=None if args.n is None else (args.n - processedReadPairsForThisLib))
                    processedReadPairsForThisLib += processedReadPairs
                    log_handle.write(
                        f"done, processed:\t{processedReadPairsForThisLib} reads\n")

                    if args.n and processedReadPairsForThisLib >= args.n:
                        break
            handle.close()
            if not args.norejects:
                rejectHandle.close()
            log_handle.write(f'Demultiplexing finished\n')
