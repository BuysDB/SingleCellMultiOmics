#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class SamTag:

    def __init__(self, tag, humanName, isPhred=False, doNotWrite=False):
        self.tag = tag
        self.humanName = humanName
        self.isPhred = isPhred
        self.doNotWrite = doNotWrite
        if len(tag) != 2:
            raise ValueError('Invalid tag length')

    def __repr__(self):
        return f'SAM TAG "{self.tag}", {self.humanName}{", is a phred score" if self.isPhred else ""} will {"not" if self.doNotWrite else ""} be written by default'


tags = [
    SamTag('SM', 'sampleName'),
    SamTag('LY', 'sequencingLibrary'),
    SamTag('EX', 'experiment'),
    SamTag('MX', 'demultiplexingStrategy'),
    SamTag('Fc', 'flowCellId'),
    SamTag('La', 'lane'),
    SamTag('Ti', 'tile'),
    SamTag('Is', 'instrument'),
    SamTag('BK', 'isBulk'),
    SamTag('XC', 'plateCoord'),
    SamTag('CX', 'clusterXpos'),
    SamTag('CY', 'clusterYpos'),
    SamTag('RP', 'readPairNumber', doNotWrite=True),
    SamTag('RN', 'runNumber'),
    SamTag('Fi', 'filterFlag'),
    SamTag('CN', 'controlNumber'),
    SamTag('aA', 'correctedIlluminaIndexSequence'),
    SamTag('aI', 'indexSequence'),
    SamTag('aa', 'rawIlluminaIndexSequence'),
    SamTag('ah', 'hammingDistanceToIndex'),
    SamTag('BC', 'barcode'),
    SamTag('CB', 'cellBarcode'),  # Same as bc
    SamTag('QX', 'barcodeQual', isPhred=True),
    SamTag('bc', 'rawBarcode'),
    SamTag('hd', 'hammingDistanceRawBcAssignedBc'),

    SamTag('bi', 'barcodeIndex'),

    #Those are insertion/deletion base qualities, where N=45, which is the default quality score for indels assigned by the BaseRecalibrator
    SamTag('BI', 'gatk_insert_deletion_base_qualities'),
    SamTag('BD', 'gatk_insert_deletion_base_qualities'),

    SamTag('QT', 'sampleBarcodeQuality', isPhred=True),
    SamTag('RX', 'umi'),
    SamTag('tu', 'second_umi'),
    SamTag('uL', 'umiLength'),
    SamTag('RQ', 'umiQual', isPhred=True),
    SamTag('BX', 'umiCorrected'),
    SamTag('BZ', 'umiCorrectedQuality', isPhred=True),
    SamTag('MI', 'molecularIdentifier'),
    SamTag('QM', 'molecularIdentifierQuality', isPhred=True),
    SamTag('DS', 'siteCoordinate'),
    SamTag('DA', 'allele'),
    SamTag('ap', 'phasedAllelicSNVs'),
    SamTag('au', 'phasednUnknownAllelicOriginSNVs'),
    SamTag('TR', 'totalRTreactions'),

    SamTag('RZ', 'recognizedSequence'),
    SamTag('RS', 'recognizedStrand'),
    SamTag('LI', 'ligationMotif'),
    SamTag('RC', 'moleculeOverseqCountIndex'),
    SamTag('af', 'associatedFragmentCount'),
    SamTag('RR', 'rejectionReason'),
    SamTag('DT', 'sourceType'),
    SamTag('EX', 'exons'),
    SamTag('IN', 'introns'),
    SamTag('GN', 'genes'),
    SamTag('JN', 'junctions'),
    SamTag('IT', 'is_transcriptome'),
    SamTag('nM', 'editDistanceToReference'),
    SamTag('NM', 'editDistanceToReference'),
    SamTag('AS', 'alignmentScore'),
    SamTag('NH', 'occurencesOfFragment'),
    SamTag('HI', 'queryHitIndex'),
    SamTag('XS', 'assignmentStatus'),
    SamTag('SQ', 'meanBaseQuality'),
    SamTag('SD', 'scarDescription'),
    SamTag('RG', 'readGroup'),
    SamTag('XA', 'alternativeAlignmentHits'),
    SamTag('MD', 'alignment'),
    SamTag('dt', 'assingedDataType'),

    SamTag('rt', 'rtReactionIndex'),
    SamTag('rd', 'rtDuplicateIndex'),

    SamTag('fS', 'fragmentSizeTag'),
    SamTag('fe', 'fragmentEndTag'),
    SamTag('fs', 'fragmentStartTag'),
    SamTag('Us', 'undigestedSiteCount'),
    SamTag('lh', 'ligationOverhangSequence'),
    SamTag('lq', 'ligationOverhangQuality', isPhred=True),

    SamTag('a1', 'Read1Clipped3primeBasesKmer'),
    SamTag('A2', 'Read2Clipped3primeBasesKmer'),
    SamTag('aQ', 'Read1Clipped3primeBasesKmerQuality', isPhred=True),
    SamTag('AQ', 'Read2Clipped3primeBasesKmerQuality', isPhred=True),

    SamTag('e1', 'Read1Clipped5primeBasesKmer'),
    SamTag('eQ', 'Read2Clipped5primeBasesKmer'),
    SamTag('eB', 'Read1Clipped5primeCycle'),
    SamTag('EB', 'Read2Clipped5primeCycle'),
    SamTag('E2', 'Read1Clipped5primeBasesKmerQuality', isPhred=True),
    SamTag('EQ', 'Read2Clipped5primeBasesKmerQuality', isPhred=True),

    SamTag('ES', 'EnzymeSequence'),  # Christoph
    SamTag('eq', 'EnzymeSequenceQualities', isPhred=True),
    SamTag('IS', 'ispcrSequence'),
    SamTag('is', 'ispcrSequenceQualities', isPhred=True),

    SamTag('H0', 'HexamerSequenceR1'),  # anna
    SamTag('H1', 'HexamerPhredScoreR1', isPhred=True),
    SamTag('H2', 'HexamerSequenceR2'),
    SamTag('H3', 'HexamerPhredScoreR2', isPhred=True),

    SamTag('XM', 'perBaseMethylationStatus'),
    SamTag('Cm', 'modified3BPContexts'),
    SamTag('Qm', 'modified4BPContexts'),
    SamTag('Cu', 'unmodified3BPContexts'),
    SamTag('Qu', 'unmodified4BPContexts'),
    SamTag('MC', 'methylatedBaseCountForFragment'),
    SamTag('uC', 'totalUnmethylated'),
    SamTag('sZ', 'Total CPG methylated'),
    SamTag('sz', 'Total CPG unmethylated'),
    SamTag('sH', 'Total CHH methylated'),
    SamTag('sh', 'Total CHH unmethylated'),
    SamTag('sX', 'Total CHG methylated'),
    SamTag('sx', 'Total CHG unmethylated'),


    SamTag('rS', 'randomPrimerSequence'),
    SamTag('rP', 'randomPrimerStart'),
    SamTag('TR', 'totalRTreactions'),
    SamTag('TF', 'totalAssociatedFragmentsIncludingOverlflow'),
    SamTag('af', 'totalAssociatedFragmentsInBam'),

    SamTag('AI', 'AssignedIntrons'),
    SamTag('AE', 'AssingedExons'),
    SamTag('iH', 'IntronHits amount of bases aligned to intron'),
    SamTag('eH', 'ExonHits amount of bases aligned to exon'),
    SamTag('SP', 'IsSpliced'),



]
