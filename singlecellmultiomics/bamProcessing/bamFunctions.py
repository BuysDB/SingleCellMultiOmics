import os
import pysam

def get_reference_from_pysam_alignmentFile(pysam_AlignmentFile, ignore_missing=False):
    try:
        for x in pysam_AlignmentFile.header.as_dict()['PG']:
            if  x.get('ID')!='bwa':
                continue
            for argument in x.get('CL').split():
                if (argument.endswith('.fa') or argument.endswith('.fasta') or argument.endswith('.fasta.gz') or argument.endswith('.fa.gz')) and (ignore_missing or os.path.exists(argument)):
                    return argument
    except Exception as e:
        pass



def sort_and_index(unsorted_path, sorted_path, remove_unsorted=False):
    """ Sort and index a bam file """
    cmd = f"""samtools sort {unsorted_path} > {sorted_path}; samtools index {sorted_path};"""
    if remove_unsorted:
        cmd = cmd +f"rm {unsorted_path};"
    os.system(cmd)


def add_readgroups_to_header( origin_bam_path, readgroups_in, target_bam_path=None ):
    """ Add the readgroups in the set readgroups to the header of origin_bam_path.

    This function first loads the header of the origin to memory.
    The supplied readgroups are added to this header.
    The new header is then exported to a SAM file. The SAM file is then
    concatenated to the original bam file.

    Args:

        origin_bam_path(str) : path to bam file to which to add readgroups to header

        readgroups_in(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
                { 'ID': ID, 'LB':library,
                'PL':platform,
                'SM':sampleLib,
                'PU':readGroup }

        target_bam_path(str) : path to write bam file including the readgrouped header to. When not supplied the output is written to the input bam file

    """

    if target_bam_path is None:
        target_bam_path = origin_bam_path

    # Create a read group dictionary
    if type(readgroups_in) is set:
        readGroupsDict={}
        for readGroup in readgroups_in:
            flowCell,lane,sampleLib = readGroup.split('.')
            try:
                library,_ = sampleLib.rsplit('_',1)
            except Exception as e:
                # the library is not part of the sample name:
                library = 'undefinedLibrary'
            readGroupsDict[readGroup] = {
                    'ID':readGroup,
                    'LB':library,
                    'PL':'ILLUMINA',
                    'SM':sampleLib,
                    'PU':readGroup}
    elif type(readgroups_in) is dict:
        readGroupsDict = readgroups_in
    else:
        raise ValueError("supply a set or dict for readgroups_in")


    with pysam.AlignmentFile(origin_bam_path, "rb") as origin:
        header = origin.header.copy()

        # Write the new header to this sam file:
        headerSamFilePath = origin_bam_path.replace('.bam','')+'.header.sam'

        hCopy = header.to_dict()
        hCopy['RG'] = list(readGroupsDict.values())

        # Write the sam file with the complete header:
        with pysam.AlignmentFile(headerSamFilePath,'w',header=hCopy):
            pass


    # Write the re-headered bam file to this path
    complete_temp_path = origin_bam_path.replace('.bam','')+'.rehead.bam'

    # Concatenate and remove origin
    rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {origin_bam_path}; }} | samtools view -b > {complete_temp_path} ;
            mv {complete_temp_path } {target_bam_path};rm {headerSamFilePath};
            """
    os.system(rehead_cmd)
