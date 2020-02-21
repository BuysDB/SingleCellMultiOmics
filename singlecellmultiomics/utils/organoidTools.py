#!/usr/env python3

import collections
import itertools
import string
import math
def cellToReplicatePlateAndPassage( cellName ):
    #' APKS3-P23-1-1s_cell18.bam'

    replicate = 'APKS'+ cellName.split('-')[0].replace('APKS','').replace('P','25')
    if 'APKS1P25' in cellName:
        passage=25
        replicate='APKS1'
        plate= '-'.join( cellName.split('_')[0].split('-')[1].replace('NLAP','').replace('L',''  ).replace('s',''))
    else:
        passage = int(cellName.split('-')[1].replace('P',''))

        plate= '-'.join( [x.replace('s',''  )  for x in cellName.split('_' )[0].split('-')[2:]])
        #plate=  '.'.join( cellName.split('_')[0].split('-')[1:] )


    return( replicate, plate, passage )


def chrom_sort_human(l):
    residual = l.replace('chr','').replace('P','').strip()
    rest = 0
    if '_' in residual:
        allele=None
        if residual.count('_')==1:
            pre, aft = residual.split('_')
        else:
            pre, allele, aft = residual.split('_')
        if allele=='A':
            rest = 0
        else:
            rest = 1
        if '-' in aft:
            rest=int(aft.split('-')[0])
        residual = pre

    try:
        return (int(residual), rest)
    except Exception as e:
        return(100, rest)

def chrom_sort_human2(l):
    residual = l.replace('chr','').replace('P','').replace('-','_').strip()
    rest = 0
    if '_' in residual:
        allele=None
        if residual.count('_')==1:
            pre, aft = residual.split('_')
        else:
            pre, allele, aft = residual.split('_')
        if allele=='A':
            rest = 0
        else:
            rest = 1
        if '-' in aft:
            rest=int(aft.split('-')[0])
        residual = pre

    try:
        return (int(residual), rest)
    except Exception as e:
        return(100, rest)


def chromIsAllelic(chrom):
    return chrom.endswith('_A') or chrom.endswith('_B')

def cellToReplicatePassagePlateAndCell( cellName ):

    if cellName.count('_')==3: # was alreaddy formatted:
        parts = cellName.split('_')
        replicate = parts[0]
        passage = int(parts[1])
        plate=parts[2]
        cell=int(parts[3].replace('cell',''))

    else:
        #' APKS3-P23-1-1s_cell18.bam'
        #for replicate in ['APKS','APK','AP','A']:
        if cellName.startswith('APKS'):
            replicate = 'APKS'+ cellName.split('-')[0].replace('APKS','')
            if 'APKS1P25' in cellName:
                replicate=replicate.replace('P','25')
        else:
            replicate = cellName.split('-')[0]

        if 'APKS1P25' in cellName:
            passage=25
            replicate='APKS1'
            plate= '-'.join( cellName.split('_')[0].split('-')[1].replace('NLAP','').replace('L',''  ).replace('s',''))

        else:
            passage = int(cellName.split('-')[1].replace('P',''))



            plate= '-'.join( [x.replace('s',''  )  for x in cellName.split('_' )[0].split('-')[2:]])
            plate = plate.replace('NLAP','').replace('L','-'  )
            #plate=  '.'.join( cellName.split('_')[0].split('-')[1:] )
        if 'cell' in cellName:
            cell = int(cellName.split('cell')[-1])
        else:
            cell = int(cellName.split('_')[-1])

    return( replicate, passage, plate, cell )


##
def bulkNameToReplicatePlateAndPassage( bulkSampleName ):
    #' APKS3-P23-1-1s_cell18.bam'

    try:
        if bulkSampleName.startswith('APKS') :
            replicate = 'APKS'+ bulkSampleName.split('-')[0].replace('APKS','').replace('P','25')
            if 'APKS1P25' in bulkSampleName:
                passage=25
                replicate='APKS1'
                plate= '-'.join( bulkSampleName.split('_')[0].split('-')[1].replace('NLAP','').replace('L',''  ).replace('s',''))
            else:
                rawPasssage = bulkSampleName.split('-')[1].replace('P','').replace('G','')
                passage=None
                try:
                    passage = int(rawPasssage)
                except Exception:
                    print(f'Could not parse {rawPasssage}')
                plate= '-'.join( [x.replace('s',''  )  for x in bulkSampleName.split('_' )[0].split('-')[2:]])
            return( replicate, plate, passage )
        else:

            replicate = None
            if bulkSampleName.startswith('APK'):
                replicate ='APK'
            elif bulkSampleName.startswith('AP'):
                replicate ='AP'
            elif bulkSampleName.startswith('A'):
                replicate ='A'
            elif bulkSampleName.startswith('WT'):
                replicate ='WT'
            else:
                raise ValueError('Unkown replicate name')

            parts = bulkSampleName.replace(replicate, '').split('_')
            passage = int(parts[-1].replace('G',''))
            plate = None
            return( replicate, plate, passage )
    except Exception as e:
        raise ValueError(f"Could not parse: {bulkSampleName}")


def cellNameToFacsName( df, library, returnPlateType=False ):
    well2index = collections.defaultdict(dict)
    index2well = collections.defaultdict(dict)
    rows = string.ascii_uppercase[:16]
    columns = range(1,25)
    for ci in range(1,385):
        i = ci-1
        rowIndex = math.floor( i/len(columns) )
        row = rows[rowIndex]
        column = columns[i%len(columns)]
        well2index[384][(f'{row}{column}',1)] = ci
        index2well[384][ci, 1] = (f'{row}{column}')

    w2i96 = {}
    for ci in range(1,97): # barcode
        i = (ci-1)*2
        initRow = math.floor( i/len(columns) )*2
        initCol = i%len(columns)
        for lib in range(1,5):
            i = (ci-1)
            column = columns[initCol+(lib==2 or lib==4)]
            row = rows[initRow+(lib>=3)]
            well2index[96][(f'{row}{column}', lib)] = ci
            index2well[96][ci, lib] = (f'{row}{column}')

    plateMap = {}
    for replicate, passage, plate, cell in df.index:

        verdict = 96 if plate.count('-')>0 else 384

        plateName = (replicate, passage, plate)
        if not plateName in plateMap:
            plateMap[plateName] = verdict
        elif plateMap[plateName]==384 and verdict==96:
            plateMap[plateName]=96
    remapped = []
    print(verdict)
    #APKS1_10_1_A1
    if returnPlateType:
        plateTypes = []
    for replicate, passage, plate, cell in df.index:
        libary =int( plate.split('-')[0])
        plateType = plateMap[(replicate, passage, plate)]
        if returnPlateType:
            plateTypes.append(plateType)
        #cell = int(cell.replace('cell',''))
        well = index2well[plateType][cell, library]
        remapped.append(f'{replicate}_{passage}_{plate}_{well}')

    if returnPlateType:
        return remapped, plateTypes

    return remapped
