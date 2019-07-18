#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import pkg_resources
import singlecellmultiomics.barcodeFileParser.barcodeFileParser as barcodeFileParser
from singlecellmultiomics.modularDemultiplexer.demultiplexingStrategyLoader import DemultiplexingStrategyLoader
import singlecellmultiomics.libraryDetection.sequencingLibraryListing as sequencingLibraryListing
import glob
import fnmatch
import os
from types import SimpleNamespace
barcode_dir = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/barcodes/')
index_dir = pkg_resources.resource_filename('singlecellmultiomics','modularDemultiplexer/indices/')

barcode_parser = barcodeFileParser.BarcodeParser(hammingDistanceExpansion=0, barcodeDirectory=barcode_dir)
index_parser =  barcodeFileParser.BarcodeParser(hammingDistanceExpansion=1, barcodeDirectory=index_dir)

dmx = DemultiplexingStrategyLoader(barcodeParser=barcode_parser,
                                   indexParser=index_parser,
                                   indexFileAlias=None)



argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""Check multiplexability of many fastq files""")

argparser.add_argument('-locations', default='.')

arguments = argparser.parse_args()
sequencing_dirs = arguments.locations.split(',')
dmxes = sorted([x.shortName  for x in dmx.demultiplexingStrategies])
print('\t'.join(['RUN','SEQ','AVO','INDEX','LIBRARY']+dmxes))


matches = []
for sdir in sequencing_dirs:
    for root, dirnames, filenames in os.walk(sdir):
        for d in dirnames:
            try:
                fp = os.path.join(root, d)
                #print('   ' + fp)
                if len(list(glob.glob(fp+'/*.fastq.gz'))) and d != 'BaseCalls':
                    matches.append(fp)
                    fastqfiles = list(glob.glob(fp+'/*.fastq.gz'))

                    args =     SimpleNamespace(replace=None,
                        fastqfiles= fastqfiles, slib=None, merge='_',
                                               dsize=1000,
                                               maxAutoDetectMethods=100,
                                               minAutoDetectPct=1
                    )
                    libraries = sequencingLibraryListing.SequencingLibraryLister(verbose=False).detect(fastqfiles,
                                                                                          args=args)
                    for library, associated_fastqs_lane in libraries.items():
                        # Obtain run id
                        run_id = '?'
                        seqid = '?'
                        index='?'
                        avo_id ='?'
                        i=None
                        dpos=None
                        for lane, reads in  associated_fastqs_lane.items():
                            parts = os.path.dirname( reads['R1'][0] ).split('/')

                            if 'Data' in parts:
                                try:
                                    i = parts.index('Data')
                                except Exception as e:
                                    i=-2

                                try:
                                    dpos = parts.index('BaseCalls')
                                except Exception as e:
                                    dpos = -1
                                    pass


                                try:
                                    run_id = parts[i-2]
                                except Exception as e:
                                    pass

                                try:
                                    seqid = parts[i-1]
                                except Exception as e:
                                    pass

                                try:
                                    index = parts[dpos+2]
                                except Exception as e:
                                    pass

                                try:
                                    avo_id = parts[dpos+1]
                                except Exception as e:
                                    pass

                                else:
                                    try:
                                        avo_id = parts[-2]
                                    except Exception as e:
                                        pass

                                    try:
                                        index = parts[-1]
                                    except Exception as e:
                                        pass

                                break

                        processedReadPairs, strategyYieldsForAllLibraries = dmx.detectLibYields(
                            {library:associated_fastqs_lane},
                            testReads=args.dsize,
                            maxAutoDetectMethods=args.maxAutoDetectMethods,
                            minAutoDetectPct=args.minAutoDetectPct)
                        print('\t'.join( [run_id,seqid,avo_id,index,library] + [ str(strategyYieldsForAllLibraries[library]['strategyYields'].get(x,0)/10) for x in dmxes]))
            except Exception as e:
                print(e)
