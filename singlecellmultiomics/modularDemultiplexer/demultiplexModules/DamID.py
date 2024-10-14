from singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods import ScatteredUmiBarcodeDemuxMethod, UmiBarcodeDemuxMethod, NonMultiplexable, IlluminaBaseDemultiplexer, phredToFastqHeaderSafeQualities
from singlecellmultiomics.modularDemultiplexer.demultiplexModules import CELSeq2_c8_u6
import re
from singlecellmultiomics.utils import reverse_complement




class DamID2(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, random_primer_read=None,random_primer_length=None, **kwargs):
        self.barcodeFileAlias = 'DamID2'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=10,
            random_primer_read=random_primer_read,
            random_primer_length=random_primer_length,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'DamID2'
        self.longName = 'DamID,  CB: 10bp UMI: 3bp'
        self.autoDetectable = True
        self.description = '48 barcodes, 3bp umi followed by 10bp barcode, the last two bases of the barcode are CA'

        self.sequenceCapture[0] = slice(
            self.umiLength + self.barcodeLength - 1,
            None)  # do map the first base of the barcode

    def demultiplex(self, records, probe=False, **kwargs):

        if probe and records[0].sequence[self.barcodeLength + self.umiLength-2:self.barcodeLength + self.umiLength] != 'CA':
            raise NonMultiplexable


        # add first 2 bases as ligation tag:
        ligation_start = self.umiLength + self.barcodeLength -2
        ligation_end = ligation_start + 2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
            self, records, probe=probe, **kwargs)

        for rec in taggedRecords:
            rec.addTagByTag(
                'lh',
                ligation_sequence,
                isPhred=False,
                make_safe=False)
            rec.addTagByTag(
                'lq',
                ligation_qualities,
                isPhred=True,
                make_safe=False)
            rec.tags['dt'] = 'DamID'

        return taggedRecords


class DamID2_NO_OVERHANG(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, random_primer_read=None,random_primer_length=None, **kwargs):
        self.barcodeFileAlias = 'DamID2_8bp'
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
        self.shortName = 'DamID2_8bp_noCA'
        self.longName = 'DamID,  CB: 8bp UMI: 3bp'
        self.autoDetectable = True
        self.description = '48 barcodes, 3bp umi followed by 8bp barcode, no CA overhang in barcodes'

        self.sequenceCapture[0] = slice(
            self.umiLength + self.barcodeLength - 1,
            None)  # do map the first base of the barcode

    def demultiplex(self, records, probe=False, **kwargs):

        #if probe and records[0].sequence[self.barcodeLength + self.umiLength-2:self.barcodeLength + self.umiLength] != 'CA':
        #    raise NonMultiplexable


        # add first 2 bases as ligation tag:
        ligation_start = self.umiLength + self.barcodeLength 
        ligation_end = ligation_start + 2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]

        taggedRecords = UmiBarcodeDemuxMethod.demultiplex(
            self, records, probe=probe, **kwargs)

        for rec in taggedRecords:
            rec.addTagByTag(
                'lh',
                ligation_sequence,
                isPhred=False,
                make_safe=False)
            rec.addTagByTag(
                'lq',
                ligation_qualities,
                isPhred=True,
                make_safe=False)
            rec.tags['dt'] = 'DamID'
        return taggedRecords
    
class DamID2_c8_u3_cs2(UmiBarcodeDemuxMethod):
    def __init__(self, barcodeFileParser, random_primer_read=None, random_primer_length=None, **kwargs):
        self.barcodeFileAlias = 'DamAndT'
        UmiBarcodeDemuxMethod.__init__(
            self,
            umiRead=0,
            umiStart=0,
            umiLength=3,
            barcodeRead=0,
            barcodeStart=3,
            barcodeLength=10,
            random_primer_read=random_primer_read,
            random_primer_length=random_primer_length,
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'DamAndT'
        self.longName = 'DamAndT,  CB: 10bp UMI: 3bp'
        self.autoDetectable = True
        self.description = 'DamID2 and CS2 Transcriptome'

        self.sequenceCapture[0] = slice(
            self.umiLength + self.barcodeLength - 2,
            None)  # do map the first two bases of the barcode
        
        self.transcriptome_demux =  CELSeq2_c8_u6(barcodeFileParser=barcodeFileParser,**kwargs)
        self.damid_demux =  DamID2(barcodeFileParser=barcodeFileParser,**kwargs)


    def demultiplex(self, records, probe=False, **kwargs):
        # Check if the supplied reads are mate-pair:
        if len(records) != 2:
            raise NonMultiplexable('Not mate pair')


        # Obtain the DAMID barcode and umi:
        try:
            damid_tagged_records = self.damid_demux.demultiplex(records, probe=probe, **kwargs)
        except NonMultiplexable as e:
            damid_tagged_records = None
            
        
        # Obtain TX barcode and umi
        try:
            tx_tagged_records = self.transcriptome_demux.demultiplex(records, probe=probe, **kwargs)
            # Prune the poly T off R1 start
            
            r1 = tx_tagged_records[0]
            
            pos = 0
            for pos, base in enumerate(r1.sequence):
                if base!='T':
                    break

            r1.sequence = r1.sequence[pos:]
            r1.qualities  = r1.qualities[pos:]

        except NonMultiplexable:
            tx_tagged_records = None

        # Check which of the demuxers yielded results
        # damid_tagged_records
        # damid_taggedRecords[0].tags.update(ud)
        # damid_taggedRecords[1].tags.update(ud)
        if tx_tagged_records is not None and damid_tagged_records is not None:
            # Both types
            tagged_records = damid_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'Ambiguous'

        elif tx_tagged_records is not None:
            tagged_records = tx_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'RNA'
        elif damid_tagged_records is not None:
            tagged_records = damid_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'DamID'
        elif tx_tagged_records is None and damid_tagged_records is None:

            raise NonMultiplexable('Not DamID or TX')
        else:
            raise ValueError('failure')

        return tagged_records


class DamID2_SCA(ScatteredUmiBarcodeDemuxMethod):
    def __init__(self, 
                 barcodeFileParser, 
                random_primer_read=None, 
                random_primer_length=None, 
                first_umi_len = 3,
                first_bc_len = 4,
                second_umi_len = 3,
                second_barcode_len = 4,
                barcode_alias=None, 
                **kwargs):
        if barcode_alias is None:
            self.barcodeFileAlias = 'DamID2_scattered_8bp'
        else:
            self.barcodeFileAlias= barcode_alias

        # eg 3-ATTG-3-GAAC[GA]
        
        self.total_mi_len = first_umi_len+first_bc_len+second_umi_len+second_barcode_len
        barcode_slices = (
            ( slice(first_umi_len,first_umi_len+first_bc_len), 
                slice(first_umi_len+first_bc_len+second_umi_len,first_umi_len+first_bc_len+second_umi_len+second_barcode_len )),  ()
        )
        umi_slices = (
            (slice(0,first_umi_len), 
                slice(first_umi_len+first_bc_len,first_umi_len+first_bc_len+second_umi_len)
                ),
            ()
        )
        capture_slices = (
            slice(self.total_mi_len,None),
            slice(None)
        )


        ScatteredUmiBarcodeDemuxMethod.__init__(
            self,
            barcode_slices = barcode_slices, # 2-len Tuple of List of slices where the cell barcode is present in the read (slice(), slice(), ...], [ slice())
            umi_slices = umi_slices, # 2-len List of slices where the UMI is present in the read [slice(), slice(), ...], [slice()] # 
            capture_slices = capture_slices, # Slices which indicate what region of read 1 and read2 to store in the final read (required)
            barcodeFileAlias=self.barcodeFileAlias,
            barcodeFileParser=barcodeFileParser,
            **kwargs)
        self.shortName = 'DamID2_3u4b3u6b'
        self.longName = 'DamID, 3bp UMI, 4bp CB, 3bp UMI, 6bp CB'
        self.autoDetectable = True
        self.description = 'DamID, starting with a 3bp UMI, 4bp CB, 3bp UMI, 6bp CB'


    def demultiplex(self, records, probe=False, **kwargs):

        #if probe and records[0].sequence[self.barcodeLength + self.umiLength-2:self.barcodeLength + self.umiLength] != 'CA':
        #    raise NonMultiplexable


        # add first 2 bases as ligation tag:
        ligation_start = self.total_mi_len
        ligation_end = self.total_mi_len+2
        ligation_sequence = records[0].sequence[ligation_start:ligation_end]
        ligation_qualities = records[0].qual[ligation_start:ligation_end]

        taggedRecords = ScatteredUmiBarcodeDemuxMethod.demultiplex(
            self, records, probe=probe, **kwargs)

        for rec in taggedRecords:
            rec.addTagByTag(
                'lh',
                ligation_sequence,
                isPhred=False,
                make_safe=False)
            rec.addTagByTag(
                'lq',
                ligation_qualities,
                isPhred=True,
                make_safe=False)

        return taggedRecords
    

class DamID2andT_SCA6(DamID2_c8_u3_cs2):
    def __init__(self, barcodeFileParser, random_primer_read=None, random_primer_length=None, **kwargs):
        
        self.transcriptome_demux =  DamID2_SCA(barcodeFileParser=barcodeFileParser,barcode_alias='CS2_scattered_8bp',
                                                first_umi_len = 3,
                                                first_bc_len = 4,
                                                second_umi_len = 3,
                                                second_barcode_len = 4,
                                               **kwargs)
        self.damid_demux =  DamID2_SCA(barcodeFileParser=barcodeFileParser,
                                       second_barcode_len=6,
                                       barcode_alias='DamID2_scattered_10bp',
                                       **kwargs)
        self.autoDetectable = True

        self.shortName = 'DamID2andT_3u4b3u6b'
        self.longName = 'DamID and transcriptome, 3bp UMI, 4bp CB, 3bp UMI, 6bp CB'
        self.autoDetectable = True
        self.description = 'DamID and transcriptome, starting with a 3bp UMI, 4bp CB, 3bp UMI, 6bp CB'
        self.indexSummary = self.transcriptome_demux.indexSummary + ', ' + self.damid_demux.indexSummary
        self.barcodeSummary = self.transcriptome_demux.barcodeSummary + ', ' + self.damid_demux.barcodeSummary

    def demultiplex(self, records, probe=False, **kwargs):
        # Check if the supplied reads are mate-pair:
        if len(records) != 2:
            raise NonMultiplexable('Not mate pair')


        # Obtain the DAMID barcode and umi:
        try:
            damid_tagged_records = self.damid_demux.demultiplex(records, probe=probe, **kwargs)
        except NonMultiplexable as e:
            damid_tagged_records = None
            
        
        # Obtain TX barcode and umi
        try:
            tx_tagged_records = self.transcriptome_demux.demultiplex(records, probe=probe, **kwargs)
            # Prune the poly T off R1 start
            
            r1 = tx_tagged_records[0]
            
            pos = 0
            for pos, base in enumerate(r1.sequence):
                if base!='T':
                    break

            r1.sequence = r1.sequence[pos:]
            r1.qualities  = r1.qualities[pos:]

        except NonMultiplexable:
            tx_tagged_records = None

        # Check which of the demuxers yielded results
        # damid_tagged_records
        # damid_taggedRecords[0].tags.update(ud)
        # damid_taggedRecords[1].tags.update(ud)
        if tx_tagged_records is not None and damid_tagged_records is not None:
            # Both types
            for tx_tagged_record, damid_tagged_record in zip(tx_tagged_records, damid_tagged_records):
                tx_tagged_record.tags.update(damid_tagged_record.tags)
            tagged_records = tx_tagged_record
        elif tx_tagged_records is not None:
            tagged_records = tx_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'RNA'
        elif damid_tagged_records is not None:
            tagged_records = damid_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'DamID'
        elif tx_tagged_records is None and damid_tagged_records is None:

            raise NonMultiplexable('Not DamID or TX')
        else:
            raise ValueError('failure')

        return tagged_records
    

class DamID2andT_SCA(DamID2_c8_u3_cs2):
    def __init__(self, barcodeFileParser, random_primer_read=None, random_primer_length=None, **kwargs):
        
        self.transcriptome_demux =  DamID2_SCA(barcodeFileParser=barcodeFileParser,barcode_alias='CS2_scattered_8bp',
                                                first_umi_len = 3,
                                                first_bc_len = 4,
                                                second_umi_len = 3,
                                                second_barcode_len = 4,
                                               **kwargs)
        
        self.damid_demux =  DamID2_SCA(barcodeFileParser=barcodeFileParser,
                                       first_umi_len = 3,
                                                first_bc_len = 4,
                                                second_umi_len = 3,
                                                second_barcode_len = 4,
                                                barcode_alias='DamID2_scattered_8bp',
                                       **kwargs
                                       )
        self.autoDetectable = True

        self.shortName = 'DamID2andT_3u4b3u4b'
        self.longName = 'DamID and transcriptome, 3bp UMI, 4bp CB, 3bp UMI, 4bp CB'
        self.autoDetectable = True
        self.description = 'DamID and transcriptome, starting with a 3bp UMI, 4bp CB, 3bp UMI, 4bp CB'
        self.indexSummary = self.transcriptome_demux.indexSummary + ', ' + self.damid_demux.indexSummary
        self.barcodeSummary = self.transcriptome_demux.barcodeSummary + ', ' + self.damid_demux.barcodeSummary

    def demultiplex(self, records, probe=False, **kwargs):
        # Check if the supplied reads are mate-pair:
        if len(records) != 2:
            raise NonMultiplexable('Not mate pair')


        # Obtain the DAMID barcode and umi:
        try:
            damid_tagged_records = self.damid_demux.demultiplex(records, probe=probe, **kwargs)
            
        except NonMultiplexable as e:
            damid_tagged_records = None
            
        
        # Obtain TX barcode and umi
        try:
            tx_tagged_records = self.transcriptome_demux.demultiplex(records, probe=probe, **kwargs)
            # Prune the poly T off R1 start
            
            r1 = tx_tagged_records[0]
            
            pos = 0
            for pos, base in enumerate(r1.sequence):
                if base!='T':
                    break

            r1.sequence = r1.sequence[pos:]
            r1.qualities  = r1.qualities[pos:]

        except NonMultiplexable:
            tx_tagged_records = None

        # Check which of the demuxers yielded results
        # damid_tagged_records
        # damid_taggedRecords[0].tags.update(ud)
        # damid_taggedRecords[1].tags.update(ud)
        if tx_tagged_records is not None and damid_tagged_records is not None:
            # Both types
            tagged_records = damid_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'Ambiguous'

        elif tx_tagged_records is not None:
            tagged_records = tx_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'RNA'
        elif damid_tagged_records is not None:
            tagged_records = damid_tagged_records
            for rec in tagged_records:
                rec.tags['dt'] = 'DamID'
        elif tx_tagged_records is None and damid_tagged_records is None:

            raise NonMultiplexable('Not DamID or TX')
        else:
            raise ValueError('failure')

        return tagged_records
    
