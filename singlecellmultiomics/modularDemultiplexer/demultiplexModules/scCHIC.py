from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import UmiBarcodeDemuxMethod, NonMultiplexable, IlluminaBaseDemultiplexer, phredToFastqHeaderSafeQualities
from singlecellmultiomics.modularDemultiplexer.demultiplexModules import CELSeq2_c8_u6
import re
from singlecellmultiomics.utils import reverse_complement
# SCCHIC using NLAIII adapter, 384 well format with 3bp UMI followed by
# "A" base


class SCCHIC_384w_c8_u3(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, random_primer_read=1,random_primer_length=6, **kwargs):
        self.barcodeFileAlias = 'maya_384NLA'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            random_primer_read=random_primer_read,
            random_primer_length=random_primer_length,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'scCHIC384C8U3'
        self.longName = 'Single cell CHIC, 384well CB: 8bp UMI: 3bp, RP: 6BP'
        self.autoDetectable = True
        self.description = '384 well format. 3bp umi followed by 8bp barcode and a single A. R2 ends with a 6bp random primer'

        self.sequenceCapture[0] = slice(
            self.barcodeLength + self.umiLength + 1,
            None)  # dont capture the first base

    def demultiplex(self, records, **kwargs):

        if kwargs.get(
                'probe') and records[0].sequence[self.barcodeLength + self.umiLength] != 'T':
            raise NonMultiplexable


        # add first 2 bases as ligation tag:
        ligation_start = self.barcodeLength + self.umiLength
        ligation_end = ligation_start + 2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
            self, records, **kwargs)

        taggedRecords[0].addTagByTag(
            'lh',
            ligation_sequence,
            isPhred=False,
            make_safe=False)
        taggedRecords[0].addTagByTag(
            'lq',
            ligation_qualities,
            isPhred=True,
            make_safe=False)
        taggedRecords[1].addTagByTag(
            'lh',
            ligation_sequence,
            isPhred=False,
            make_safe=False)
        taggedRecords[1].addTagByTag(
            'lq',
            ligation_qualities,
            isPhred=True,
            make_safe=False)
        #taggedRecords[0].sequence = taggedRecords[0].sequence[1:]
        #taggedRecords[0].qualities = taggedRecords[0].qualities[1:]
        return taggedRecords

class SCCHIC_384w_c8_u3_direct_ligation(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, random_primer_read=None,random_primer_length=None, **kwargs):
        self.barcodeFileAlias = 'maya_384NLA'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            random_primer_read=random_primer_read,
            random_primer_length=random_primer_length,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'scCHIC384C8U3l'
        self.longName = 'Single cell CHIC, 384well CB: 8bp UMI: 3bp, no RP'
        self.autoDetectable = False
        self.description = '384 well format. 3bp umi followed by 8bp barcode and a single A. R2 does not contain a random primer'

        self.sequenceCapture[0] = slice(
            self.barcodeLength + self.umiLength + 1,
            None)  # dont capture the first base

    def demultiplex(self, records, **kwargs):

        if kwargs.get(
                'probe') and records[0].sequence[self.barcodeLength + self.umiLength] != 'T':
            raise NonMultiplexable


        # add first 2 bases as ligation tag:
        ligation_start = self.barcodeLength + self.umiLength
        ligation_end = ligation_start + 2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
            self, records, **kwargs)

        taggedRecords[0].addTagByTag(
            'lh',
            ligation_sequence,
            isPhred=False,
            make_safe=False)
        taggedRecords[0].addTagByTag(
            'lq',
            ligation_qualities,
            isPhred=True,
            make_safe=False)
        taggedRecords[1].addTagByTag(
            'lh',
            ligation_sequence,
            isPhred=False,
            make_safe=False)
        taggedRecords[1].addTagByTag(
            'lq',
            ligation_qualities,
            isPhred=True,
            make_safe=False)
        #taggedRecords[0].sequence = taggedRecords[0].sequence[1:]
        #taggedRecords[0].qualities = taggedRecords[0].qualities[1:]
        return taggedRecords


class SCCHIC_384w_c8_u3_direct_ligation_SINGLE_END(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, random_primer_read=None,random_primer_length=None, **kwargs):
        self.barcodeFileAlias = 'maya_384NLA'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            random_primer_read=None,
            random_primer_length=None,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'scCHIC384C8U3se'
        self.longName = 'Single cell CHIC, 384well CB: 8bp UMI: 3bp, single end, no RP'
        self.autoDetectable = True
        self.description = '384 well format. 3bp umi followed by 8bp barcode and a single A. No read 2'

        self.sequenceCapture[0] = slice(
            self.barcodeLength + self.umiLength + 1,
            None)  # dont capture the first base

    def demultiplex(self, records, **kwargs):

        if kwargs.get(
                'probe') and records[0].sequence[self.barcodeLength + self.umiLength] != 'T':
            raise NonMultiplexable

        if len(records) != 1:
            raise NonMultiplexable

        # add first 2 bases as ligation tag:
        ligation_start = self.barcodeLength + self.umiLength
        ligation_end = ligation_start + 2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
            self, records, **kwargs)

        taggedRecords[0].addTagByTag(
            'lh',
            ligation_sequence,
            isPhred=False,
            make_safe=False)
        taggedRecords[0].addTagByTag(
            'lq',
            ligation_qualities,
            isPhred=True,
            make_safe=False)

        return taggedRecords

class SCCHIC_384w_c8_u3_pdt(IlluminaBaseDemultiplexer):

    def __init__(
            self,
            barcodeFileParser=None,
            indexFileParser=None,
            indexFileAlias='illumina_merged_ThruPlex48S_RP',

            **kwargs):

        IlluminaBaseDemultiplexer.__init__(
            self,
            indexFileParser=indexFileParser,
            indexFileAlias=indexFileAlias)

        self.description = '384 well format, mixed transcriptome and CHiC. scCHiC: 3bp umi followed by 8bp barcode and a single A. R2 ends with a 6bp random primer. Transcriptome: cs2 + template switching oligo'
        self.shortName = 'CHICTV'

        self.autoDetectable = False

        # The demultiplexer used for the chic reads:
        self.chic_demux =  SCCHIC_384w_c8_u3(barcodeFileParser=barcodeFileParser,random_primer_read=None,**kwargs)

        self.barcodeSummary = self.chic_demux.barcodeSummary
        self.longName = f'{self.chic_demux.longName} and TV primer '

    def __repr__(self):
        return f'{self.longName} {self.description}'

    def demultiplex(self, records, **kwargs):

        # Check if the supplied reads are mate-pair:
        if len(records) != 2:
            raise NonMultiplexable('Not mate pair')

        # Check if the TSO oligo is present..
        if 'AGACTCTTT' in records[0].sequence:


            taggedRecords = self.chic_demux.demultiplex(records,**kwargs)


            if not 'AGACTCTTT' in taggedRecords[0].sequence:
                raise NonMultiplexable('No match to transcriptome or CHiC')

            oligo_position = taggedRecords[0].sequence.index('AGACTCTTT')
            umi_start = max(0,oligo_position-6)
            umi = taggedRecords[0].sequence[umi_start:oligo_position]

            # Set umi :
            for r in taggedRecords:
                r.tags['tu'] = umi
                r.tags['MX'] = "CTV"

            # Clip the read down:
            taggedRecords[0].sequence = taggedRecords[0].sequence[:oligo_position]
            taggedRecords[0].qualities = taggedRecords[0].qualities[:oligo_position]

            return taggedRecords
        raise NonMultiplexable('No match to transcriptome or CHiC')

class SCCHIC_384w_c8_u3_cs2(UmiBarcodeDemuxMethod):

    def __init__(self, barcodeFileParser, random_primer_read=None,random_primer_length=None, **kwargs):
        self.barcodeFileAlias = 'maya_384NLA'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=8,
            random_primer_read=random_primer_read,
            random_primer_length=random_primer_length,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)

        self.description = '384 well format, mixed transcriptome and CHiC. scCHiC: 3bp umi followed by 8bp barcode and a single A. R2 has no random primer. Transcriptome is VASA'
        self.shortName = 'TCHIC'

        self.autoDetectable = False

        # The demultiplexer used for the transcriptome reads:
        self.transcriptome_demux =  CELSeq2_c8_u6(barcodeFileParser=barcodeFileParser,**kwargs)

        # Contains expected bleedthrough sequence
        self.id_to_cs2_barcode = { v:k + 'TTTTT' for k,v in barcodeFileParser['celseq2'].items() }
        assert len(self.id_to_cs2_barcode)>0
        self.tx_umi_len = 6

        # The demultiplexer used for the chic reads:
        self.chic_demux =  SCCHIC_384w_c8_u3_direct_ligation(barcodeFileParser=barcodeFileParser,**kwargs)

        self.barcodeSummary = f'{self.chic_demux.barcodeSummary} and {self.transcriptome_demux.barcodeSummary}'
        self.longName = f'{self.chic_demux.longName} and {self.transcriptome_demux.longName}'

        self.poly_length=10
        self.poly_A = self.poly_length*'A'
        self.poly_T = self.poly_length*'T'
        self.poly_G = self.poly_length*'G'

        self.r2_trimmer = re.compile('[GA]*$')

        self.sequenceCapture[0] = slice(
            self.barcodeLength + self.umiLength + 1,
            None)  # dont capture the first base


    def trim_r2(self, sequence, qualities ):

        start = sequence.find(self.poly_A)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        start = sequence.find(self.poly_G)
        if start != -1:
            sequence, qualities = sequence[:start], qualities[:start]

        # Trim any trailing A and G bases from the end and # Trim down 3 bases
        sequence = self.r2_trimmer.sub('',sequence)[:-3]
        qualities = qualities[:len(sequence)]


        return sequence, qualities


    def __repr__(self):
        return f'{self.longName} {self.description}'

    def extract_vasa_umi(self, sequence, vasa_barcode):
        # Extract vasa umi from sequence given the vasa barcode
        # Returns None when the it cannot be extracted
        vasa_umi_end = sequence.find(vasa_barcode)
        vasa_umi_start = max(0,vasa_umi_end - self.tx_umi_len)
        if vasa_umi_end - vasa_umi_start > 0:
            return sequence[vasa_umi_start:vasa_umi_end]
        return None

    def demultiplex(self, records, **kwargs):

        # Check if the supplied reads are mate-pair:
        if len(records) != 2:
            raise NonMultiplexable('Not mate pair')


        # add first 2 bases as ligation tag:
        ligation_start = self.barcodeLength + self.umiLength
        ligation_end = ligation_start + 2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]
        # Obtain the chic barcode and umi:
        try:
            taggedRecords = UmiBarcodeDemuxMethod.demultiplex(self,records, **kwargs)
        except NonMultiplexable:
            raise

        # Trim ligation motif:
        # add first 2 bases as ligation tag:
        ud = {
            'lh':ligation_sequence,
            'lq':phredToFastqHeaderSafeQualities(ligation_qualities),
            'dt':'CHIC'
        }

        taggedRecords[0].tags.update(ud)
        taggedRecords[1].tags.update(ud)

        # Check for tx contamination:
        expected_barcode = self.id_to_cs2_barcode.get( taggedRecords[0].tags['bi'] )
        if expected_barcode is None:
            raise ValueError(taggedRecords[0].tags['bi'])
        if expected_barcode in taggedRecords[0].sequence or expected_barcode in reverse_complement(taggedRecords[1].sequence):
            # tx contaminant:
            # Obtain vasa UMI:
            if expected_barcode in taggedRecords[0].sequence:
                # NNNNNNNNN [UMI] BC POLY T
                tx_umi = self.extract_vasa_umi(taggedRecords[0].sequence, expected_barcode)
            else:
                rc2 = reverse_complement(taggedRecords[1].sequence)
                tx_umi = self.extract_vasa_umi(rc2, expected_barcode)

            # Trim read2 down
            taggedRecords[1].sequence, taggedRecords[1].qualities = self.trim_r2(taggedRecords[1].sequence, taggedRecords[1].qualities)
            for r in taggedRecords:
                r.tags['dt'] = 'VASA'
                if tx_umi is not None:
                    r.tags['rx'] = tx_umi

                #r.tags['MX'] = 'CS2'
            return taggedRecords
        elif 'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[0].sequence or  'TTTTTTTTTTTTTTTTTTTTTTT' in taggedRecords[1].sequence:
            raise NonMultiplexable('PolyT')
        elif 'AGTCCGACGAT' in taggedRecords[0].sequence[:30] or 'GTTCTACAGT' in taggedRecords[0].sequence[:30] or 'TAATACGACTCACTATAGGG' in taggedRecords[0].sequence:

            # desc = 'none?'
            # if 'AGTCCGACGAT' in taggedRecords[0].sequence[:30]:
            #     desc = 'AGTCCGACGAT'
            # elif 'GTTCTACAGT' in taggedRecords[0].sequence[:30]:
            #     desc = 'GTTCTACAGT'
            # elif 'TAATACGACTCACTATAGGG' in taggedRecords[0].sequence:
            #     desc = 'TAATACGACTCACTATAGGG'

            for r in taggedRecords:
                r.tags['dt'] = 'VASA'
                r.tags['RR'] = 'T7_found'
                #r.tags['TR'] = desc

        return taggedRecords
