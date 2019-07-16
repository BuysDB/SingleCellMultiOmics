#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import collections
import glob
import sys
from colorama import Fore
from colorama import Back
from colorama import Style
import importlib
import inspect
import traceback
import singlecellmultiomics.modularDemultiplexer.demultiplexModules as dm
import singlecellmultiomics.fastqProcessing.fastqIterator as fastqIterator
from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import NonMultiplexable,IlluminaBaseDemultiplexer


class DemultiplexingStrategyLoader:
    def __init__(self, barcodeParser, moduleSearchDir= 'demultiplexModules', indexParser=None, ignoreMethods=None,indexFileAlias=None):
        package = f'singlecellmultiomics.modularDemultiplexer.{moduleSearchDir}'
        moduleSearchPath = os.path.join( os.path.dirname(os.path.realpath(__file__)), moduleSearchDir).replace('/./','/')
        self.barcodeParser = barcodeParser
        self.indexParser = indexParser
        moduleSearchPath = moduleSearchPath
        print(f'{Style.DIM}Current script location: {__file__}')
        #print(f'Searchdir: {moduleSearchDir}')
        #print(f'Looking for modules in {moduleSearchPath}{Style.RESET_ALL}')
        self.demultiplexingStrategies = []
        self.demux_classes = [
            dm.CELSeq2_c8_u8,
            dm.MSPJI_c8_u3,
            dm.ScartraceR1,
            dm.CELSeq2_c8_u8_NNLAIII,
            dm.ScartraceR2,
            dm.NLAIII_384w_c8_u3,
            dm.CELSeq1_c8_u4,
            dm.NLAIII_96w_c8_u3,
            dm.Nla_384w_u8_c8_ad3_is15,
            dm.CELSeq2_c8_u6,
            IlluminaBaseDemultiplexer,
            dm.CELSeq2_c8_u6_swapped_reads,
            dm.SCCHIC_384w_c8_u3
        ]
        for c in self.demux_classes:
            self.demultiplexingStrategies.append(c(
            barcodeFileParser=barcodeParser,
            indexFileParser=indexParser,
            indexFileAlias=indexFileAlias)
            )

        """
        for modulePath in glob.glob(f'{moduleSearchPath}/*.py'):
            #print(f"Found {modulePath}")
            try:
                module = (modulePath.split('/')[-1].replace('.py',''))
                #print(f"Module:{module}, package:{package}")

                if ignoreMethods is not None and module in ignoreMethods:
                    print(f"{Style.DIM}Ignoring demultiplex method {module}, use -ignoreMethods none to re-enable{Style.RESET_ALL}")
                    continue
                if module=='__init__':
                    continue
 #modulePath.replace('\\','/').replace('/','.').replace('..','').replace('.py','').lstrip('.').split('.')[-1]

                loadedModule = importlib.import_module(f'.{module}', package)
                # Only obtain classes defined in the module, not imported ones:
                is_class_member = lambda member: inspect.isclass(member) and member.__module__ == f'{package}.{module}'
                for className, classDetails in inspect.getmembers(sys.modules[f'{package}.{module}'], is_class_member):
                    # Obtain a handle to the class and instatiate the strategy
                    if 'Base_' in className:
                        continue

                    class_ = getattr(loadedModule, className)
                    initiatedDemultiplexingStrategy = class_( barcodeFileParser=barcodeParser, indexFileParser=indexParser,indexFileAlias=indexFileAlias)
                    self.demultiplexingStrategies.append(initiatedDemultiplexingStrategy)
                    #print(initiatedDemultiplexingStrategy.name)

            except Exception as e:

                print(f"{Fore.RED}{Style.BRIGHT}FAILED LOADING {module} at {modulePath}\nException: {e}{Style.RESET_ALL}\nTraceback for the error:\n")
                import traceback
                traceback.print_exc()

                from os import stat
                from pwd import getpwuid

                print(f'Contact {Style.BRIGHT}%s{Style.RESET_ALL} for help\n' % getpwuid(stat(modulePath).st_uid).pw_name)
                print('The error only affects this module.\nProceeding to load more modules...\n')
        """
    def getSelectedStrategiesFromStringList(self, strList):
        selectedStrategies = []

        resolved = {part:False for part in strList}
        for strategy in self.demultiplexingStrategies:
            if strategy.shortName in strList:
                selectedStrategies.append(strategy)
                print('Selected strategy %s' % strategy)
                resolved[strategy.shortName] = True


        if any( [v is False for v in resolved.values()]):
            for strat in strList:
                if resolved[strat] is False:
                    print(f'{Fore.RED}Could not resolve {strat}{Style.RESET_ALL}')
                    print('Available:')
                    for s in self.demultiplexingStrategies:
                        print(s.shortName)
                    raise ValueError(f'Strategy {strat} not found')

            raise ValueError()
        return selectedStrategies

    def list(self):
        print(f"{Style.BRIGHT}Available demultiplexing strategies:{Style.RESET_ALL}")
        #print('Name, alias, will be auto detected, description')
        for strategy in self.demultiplexingStrategies:

            try:
                print(f'{Style.BRIGHT}{strategy.shortName}{Style.RESET_ALL}\t{strategy.longName}\t' + (f'{Fore.GREEN}Will be autodetected' if strategy.autoDetectable else f'{Fore.RED}Will not be autodetected')+Style.RESET_ALL + Style.DIM + f' {strategy.barcodeFileParser.getTargetCount(strategy.barcodeFileAlias) if hasattr(strategy,"barcodeFileParser") else "NA"} targets\n '+ Style.DIM + strategy.description +'\n'+strategy.getParserSummary() + Style.RESET_ALL +'\n' )
            except Exception as e:
                print(f"{Fore.RED}{Style.BRIGHT}Error in: {strategy.shortName}\nException: {e}{Style.RESET_ALL}\nTraceback for the error:\n")
                import traceback
                traceback.print_exc()
                from os import stat
                from pwd import getpwuid
                try:
                    modulePath = sys.modules[strategy.__module__].__file__

                    print(f'Contact {Style.BRIGHT}%s{Style.RESET_ALL} for help\n' % getpwuid(stat(modulePath).st_uid).pw_name)
                    print('The error only affects this module.\nProceeding to load more modules...\n')
                except Exception as e:
                    pass

    def getAutodetectStrategies(self):
        return [strategy for strategy in self.demultiplexingStrategies if  strategy.autoDetectable ]

    def getDemultiplexingSelectedStrategies(self):
        if self.selectedStrategies is None:
            raise ValueError('No strategies selected')
        return self.selectedStrategies

    def demultiplex(self, fastqfiles, maxReadPairs=None, strategies=None, library=None, targetFile=None, rejectHandle=None,probe=None):

        useStrategies = strategies if strategies is not None else self.getAutodetectStrategies()
        strategyYields = collections.Counter()
        processedReadPairs=0
        baseDemux = IlluminaBaseDemultiplexer(indexFileParser=self.indexParser, barcodeParser=self.barcodeParser,probe=probe)

        for processedReadPairs, reads in enumerate(fastqIterator.FastqIterator(*fastqfiles)):
            for strategy in useStrategies:
                try:
                    recodedRecords = strategy.demultiplex(reads, library=library,probe=probe)

                    if targetFile is not None:
                        targetFile.write( recodedRecords )

                except NonMultiplexable:
                    #print('NonMultiplexable')
                    if rejectHandle is not None:
                        try:
                            rejectHandle.write( baseDemux.demultiplex(reads, library=library) )
                        except NonMultiplexable as e:
                            rejectHandle.write( reads )

                    continue
                except Exception as e:
                    print( traceback.format_exc() )
                    print(f'{Fore.RED}Fatal error. While demultiplexing strategy {strategy.longName} yielded an error, the error message was: {e}')
                    print('The read(s) causing the error looked like this:')
                    for read in reads:
                        print(str(read))
                    print(Style.RESET_ALL)
                #print(recodedRecord)
                strategyYields[strategy.shortName]+=1
            if ( maxReadPairs is not None and (1+processedReadPairs)>=maxReadPairs):
                break
        return processedReadPairs+1,strategyYields



    def detectLibYields(self, libraries, strategies=None, testReads=100000,maxAutoDetectMethods=1,minAutoDetectPct=5):
        libYields = dict()

        for lib, lanes in libraries.items():
            for lane, readPairs in lanes.items():

                for readPair in readPairs:
                    if len(readPairs)==1:
                        processedReadPairs, strategyYields = self.demultiplex( [readPairs['R1'][0]],maxReadPairs=testReads,strategies=strategies , probe=True )
                    elif len(readPairs)==2:
                        processedReadPairs, strategyYields  =  self.demultiplex( (readPairs['R1'][0], readPairs['R2'][0] ),maxReadPairs=testReads,strategies=strategies , probe=True  )
                    else:
                        raise ValueError('Error: %s' % readPairs.keys())

                print(f'Report for {lib}:')
                self.strategyYieldsToFormattedReport( processedReadPairs, strategyYields,maxAutoDetectMethods=maxAutoDetectMethods,minAutoDetectPct=minAutoDetectPct)
                libYields[lib]= {'processedReadPairs':processedReadPairs, 'strategyYields':strategyYields }
                break
        return processedReadPairs, libYields

    def strategyYieldsToFormattedReport(self, processedReadPairs, strategyYields, selectedStrategies=None,maxAutoDetectMethods=1,minAutoDetectPct=5):
        print(f'Analysed {Style.BRIGHT}{processedReadPairs}{Style.RESET_ALL} read pairs')

        if selectedStrategies is None:
            selectedStrategies = {}
        #selectedStrategies = self.selectedStrategiesBasedOnYield(processedReadPairs, strategyYields)
        for i,(strategy, strategyYield) in enumerate(strategyYields.most_common()):
            yieldRatio = strategyYield/(0.001+processedReadPairs)
            print( ( Style.BRIGHT+Fore.GREEN if ((strategy in selectedStrategies) or i<maxAutoDetectMethods) else (Fore.YELLOW if yieldRatio*100>=minAutoDetectPct else Style.DIM)) +  f'\t {strategy}:%.2f%%{Style.RESET_ALL}'% (100.0*yieldRatio))

    def selectedStrategiesBasedOnYield(self, processedReadPairs, strategyYields, maxAutoDetectMethods=1, minAutoDetectPct=0.05):
        selectedStrategies = []
        for strategy, strategyYield in strategyYields.most_common(maxAutoDetectMethods):
            yieldRatio = strategyYield/(0.001+processedReadPairs)*100.0
            if yieldRatio>=minAutoDetectPct:
                selectedStrategies.append(strategy)
        return selectedStrategies
