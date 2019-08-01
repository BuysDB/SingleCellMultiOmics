import os
def get_reference_from_pysam_alignmentFile(pysam_AlignmentFile, ignore_missing=False):
    try:
        for x in pysam_AlignmentFile.header.as_dict()['PG']:
            if  x.get('ID')!='bwa':
                continue
            for argument in x.get('CL').split():
                if (argument.endswith('.fa') or argument.endswith('.fasta') or argument.endswith('.fasta.gz') or argument.endswith('.fa.gz')) and (ignore_missing or os.path.exist(argument)):
                    return argument
    except Exception as e:
        pass
