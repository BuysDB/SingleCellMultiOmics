import os
import pysam
import time
import contextlib
from shutil import which
from singlecellmultiomics.utils import BlockZip
import uuid
import os
from collections import defaultdict



def verify_and_fix_bam(bam_path):
    """
    Check if the bam file is not truncated and indexed.
    If not, apply index

    Args:
        bam_path(str) : path to bam file

    Raises:
        ValueError : when the file is corrupted and a fix could not be applied
    """
    index_missing = False
    if not os.path.exists(bam_path):
        raise ValueError(f"The bam file {bam_path} does not exist")

    with pysam.AlignmentFile(bam_path, "rb") as alignments:
        if alignments.check_truncation():
            raise ValueError(f"The bam file {bam_path} is truncated")

        # This raises and error:
        try:
            if not alignments.check_index():
                # Try to index the input file..
                print(
                    f"The bam file {bam_path} does not have an index, attempting an index build ..")
                index_missing = True
        except Exception as e:
            index_missing = True

        if index_missing:
            pysam.index(bam_path)

    if index_missing:
        with pysam.AlignmentFile(bam_path, "rb") as alignments:
            if not alignments.check_index():
                raise ValueError(
                    f'The file {bam_path} is not sorted or damaged in some way')


def _get_samples_from_bam(handle):
    """Get a list of samples present in the bam_file
    private: please use get_samples_from_bam()

    Args:
        bam_file(pysam.AlignmentFile) : path to bam file or pysam object

    Returns:
        samples (set) : set containing all sample names

    """
    return set([entry['SM'] for entry in handle.header.as_dict()['RG']])

def get_sample_to_read_group_dict(bam):
    """ Obtain a dictionary containing {'sample name' : ['read groupA', 'read group B'], ...}
        Args:
            bam_file(pysam.AlignmentFile) or path to bam file or pysam object
    """
    if isinstance(bam, str):
        with pysam.AlignmentFile(bam) as pysam_AlignmentFile_handle:
            return _get_sample_to_read_group_dict(pysam_AlignmentFile_handle)
    elif isinstance(bam, pysam.AlignmentFile):
        return _get_sample_to_read_group_dict(bam)

    else:
        raise ValueError(
            'Supply either a path to a bam file or pysam.AlignmentFile object')


def _get_sample_to_read_group_dict(handle):
    """ Obtain a dictionary containing {'sample name' : ['read groupA', 'read group B']}
        Args:
            bam_file(pysam.AlignmentFile) : path to bam file or pysam object
    """

    sample_to_read_group_dict = defaultdict(list)
    # Traverse the read groups;
    for entry in handle.header.as_dict()['RG']:
        sample_to_read_group_dict[ entry['SM'] ].append(entry['ID'])
    return sample_to_read_group_dict


def get_samples_from_bam(bam):
    """Get a list of samples present in the bam_file

    Args:
        bam_file(str) or pysam.AlignmentFile : path to bam file or pysam object

    Returns:
        samples (set) : set containing all sample names

    """

    if isinstance(bam, str):
        with pysam.AlignmentFile(bam) as pysam_AlignmentFile_handle:
            return _get_samples_from_bam(pysam_AlignmentFile_handle)
    elif isinstance(bam, pysam.AlignmentFile):
        return _get_samples_from_bam(bam)

    else:
        raise ValueError(
            'Supply either a path to a bam file or pysam.AlignmentFile object')



def get_reference_from_pysam_alignmentFile(
        pysam_AlignmentFile, ignore_missing=False):
    """Extract path to reference from pysam handle

    Args:
        pysam_AlignmentFile (pysam.AlignmentFile)
        ignore_missing(bool) : Check if the file exists, if not return None
    Returns:
        path : path to bam file (if exists or ignore_missing is supplied) or None
    """
    try:
        for x in pysam_AlignmentFile.header.as_dict()['PG']:
            if x.get('ID') != 'bwa':
                continue
            for argument in x.get('CL').split():
                if (argument.endswith('.fa') or argument.endswith('.fasta') or argument.endswith(
                        '.fasta.gz') or argument.endswith('.fa.gz')) and (ignore_missing or os.path.exists(argument)):
                    return argument
    except Exception as e:
        pass


def get_reference_path_from_bam(bam, ignore_missing=False):
    """Extract path to reference from bam file

    Args:
        bam(str) or pysam.AlignmentFile : path to bam file
        ignore_missing(bool) : Check if the file exists, if not return None
    Returns:
        path : path to bam file (if exists or ignore_missing is supplied) or None
    """
    if isinstance(bam, str):
        with pysam.AlignmentFile(bam) as pysam_AlignmentFile_handle:
            return get_reference_from_pysam_alignmentFile(
                pysam_AlignmentFile_handle, ignore_missing=ignore_missing)
    elif isinstance(bam, pysam.AlignmentFile):
        return get_reference_from_pysam_alignmentFile(
            bam, ignore_missing=ignore_missing)

    else:
        raise ValueError(
            'Supply either a path to a bam file or pysam.AlignmentFile object')


@contextlib.contextmanager
def sorted_bam_file(
        write_path,
        origin_bam=None,
        header=None,
        read_groups=None,
        local_temp_sort=True):
    """ Get writing handle to a sorted bam file

    Args:
        write_path (str) : write to a  bam file at this path

        origin_bam (pysam.AlignmentFile ) : bam file to copy
         header for or
        header (dict) : header for the bam file to write

        read_groups(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
                { 'ID': ID, 'LB':library,
                'PL':platform,
                'SM':sampleLib,
                'PU':readGroup }

        local_temp_sort(bool) : create temporary files in current directory

    Example:
        >>> # This example assumes molecules are generated by `molecule_iterator`
        >>> read_groups = set() # Store unique read groups in this set
        >>> with sorted_bam_file('test_output.bam', header=input_header,read_groups=read_groups) as out:
        >>>     for molecule in molecule_iterator:
        >>>         molecule.write_tags()
        >>>         molecule.write_pysam(out)
        >>>         for fragment in molecule:
        >>>             read_groups.add(fragment.get_read_group())
        test_output.bam will be written, with read groups defined, sorted and indexed.


    Write some pysam reads to a sorted bam file
    Example:
        >>> import pysam
        >>> from singlecellmultiomics.bamProcessing import sorted_bam_file
        >>> test_sam = pysam.AlignmentFile('test.sam','w',reference_names=['chr1','chr2'],reference_lengths=[1000,1000])
        >>> read_A = pysam.AlignedSegment(test_sam.header)
        >>> read_A.reference_name = 'chr2'
        >>> read_A.reference_start = 100
        >>> read_A.query_sequence = 'TTGCA'
        >>> read_A.query_name= 'READ_A'
        >>> read_A.cigarstring = '5M'
        >>> read_A.query_qualities = [30] * len(read_A.query_sequence)
        >>> read_A.set_tag('RG','HVKCCBGXB.4.MYLIBRARY_1')
        >>> read_B = pysam.AlignedSegment(test_sam.header)
        >>> read_B.reference_name = 'chr1'
        >>> read_B.reference_start = 100
        >>> read_B.query_sequence = 'ATCGGG'
        >>> read_B.cigarstring = '6M'
        >>> read_B.query_name= 'READ_B'
        >>> read_B.query_qualities = [30] * len(read_B.query_sequence)
        >>> read_B.set_tag('RG','HVKCCBGXB.4.MYLIBRARY_2')
        >>> read_groups = set(( 'HVKCCBGXB.4.MYLIBRARY_2','HVKCCBGXB.4.MYLIBRARY_1'))
        >>> with sorted_bam_file('out.bam', header=test_sam.header,read_groups=read_groups) as out:
        >>>     out.write(read_A)
        >>>     out.write(read_B)
        Results in the bam file:
        @HD	VN:1.6	SO:coordinate
        @SQ	SN:chr1	LN:1000
        @SQ	SN:chr2	LN:1000
        @RG	ID:HVKCCBGXB.4.MYLIBRARY_2	SM:MYLIBRARY_2	LB:MYLIBRARY	PU:HVKCCBGXB.4.MYLIBRARY_2	PL:ILLUMINA
        @RG	ID:HVKCCBGXB.4.MYLIBRARY_1	SM:MYLIBRARY_1	LB:MYLIBRARY	PU:HVKCCBGXB.4.MYLIBRARY_1	PL:ILLUMINA
        READ_B	0	chr1	101	0	6M	*	0	0	ATCGGG	??????	RG:Z:HVKCCBGXB.4.MYLIBRARY_2
        READ_A	0	chr2	101	0	5M	*	0	0	TTGCA	?????	RG:Z:HVKCCBGXB.4.MYLIBRARY_1


    """
    unsorted_path = None
    unsorted_alignments = None
    target_dir = None

    unsorted_path = f'{write_path}.unsorted'
    # Create output folder if it does not exists
    target_dir = os.path.dirname(unsorted_path)
    if not os.path.exists(target_dir) and len(target_dir)>0 and target_dir!='.':
        try:
            os.makedirs(target_dir, exist_ok=True)
        except Exception as e:
            pass

    if header is not None:
        pass
    elif origin_bam is not None:
        header = origin_bam.header.copy()
    else:
        raise ValueError("Supply a header or origin_bam object")


    try: # Try to open the temp unsorted file:
        unsorted_alignments = pysam.AlignmentFile(
            unsorted_path, "wb", header=header)
    except Exception:
        # Raise when this fails
        raise

    # Yield a handle to the alignments,
    # this handle will be released when the handle runs out of scope
    yield unsorted_alignments
    unsorted_alignments.close()

    if read_groups is not None:
        add_readgroups_to_header(unsorted_path, read_groups)

    # Write, sort and index
    sort_and_index(
        unsorted_path,
        write_path,
        remove_unsorted=True,
        local_temp_sort=local_temp_sort)


def write_program_tag(input_header,
                      program_name,
                      command_line,
                      description,
                      version
                      ):
    """Write Program Tag to bam file header
    Args:
        input_header  (dict): header to write PG tag to

        program_name (str) : value to write to PN tag

        command_line (str) : value to write to CL tag

        version (str) : value to write to VN tag

        description (str) : value to write to DS tag

    """
    if 'PG' not in input_header:
        input_header['PG'] = []

    input_header['PG'].append({
        'PN': program_name,
        'CL': command_line,
        'VN': version,
        'DS': description
    })


def sort_and_index(
        unsorted_path,
        sorted_path,
        remove_unsorted=False,
        local_temp_sort=True):
    """ Sort and index a bam file
    Args:
        unsorted_path (str) : path to unsorted bam file

        sorted_path (str) : write sorted file here

        remove_unsorted (bool) : remove the unsorted file

        local_temp_sort(bool): create temporary files in target directory
    Raises:
        SamtoolsError when sorting or indexing fails
    """
    if local_temp_sort:
        base_directory = os.path.abspath( os.path.dirname(sorted_path) )
        prefix = f'TMP.{uuid.uuid4()}'
        temp_path_first = f'{base_directory}/{prefix}'
        if temp_path_first.startswith('/TMP'):
            # Perform sort in current directory
            temp_path_first = f'./{prefix}'


        # Try to sort at multiple locations, if sorting fails try the next until all locations have been tried
        temp_paths =  [temp_path_first, f'/tmp/{prefix}', f'./{prefix}' ]
        for i,temp_path in enumerate(temp_paths):
            failed=False
            try:
                pysam.sort(
                    '-o',
                    sorted_path,
                    '-T',
                    f'{temp_path}', # Current directory with a random prefix
                    unsorted_path,
                )
            except Exception as e:
                failed = True
                if i==len(temp_paths)-1:
                    raise

            if not failed:
                break
    else:
        pysam.sort("-o", sorted_path, unsorted_path)
    pysam.index(sorted_path)
    if remove_unsorted:
        os.remove(unsorted_path)


class MapabilityReader():

    def __init__(self, mapability_safe_file_path):
        self.mapability_safe_file_path = mapability_safe_file_path
        self.handle = BlockZip(mapability_safe_file_path, 'r')

        # Todo: exit statements

    def site_is_mapable(self, contig, ds, strand):
        """ Obtain if a restriction site is mapable or not
        Args:
            contig (str) : contig of site to look up
            ds (int) : zero based coordinate of site to look up
            strand (bool) : strand of site to look up (False: FWD, True: REV)

        Returns:
            site_is_mapable (bool) : True when the site is uniquely mapable, False otherwise
        """
        if self.handle[contig, ds, strand] == 'ok':
            return True
        return False


def GATK_indel_realign(origin_bam, target_bam,
                       contig, region_start, region_end,
                       known_variants_vcf_path,
                       # realignerTargetCreatorArgs=None,
                       # indelRealignerArgs=None,
                       gatk_path='GenomeAnalysisTK.jar',
                       interval_path=None,
                       java_cmd='java -jar -Xmx40G -Djava.io.tmpdir=./gatk_tmp',
                       reference=None,
                       interval_write_path=None
                       ):
    """
    Re-align a specified region in a bam file using GenomeAnalysisTK

    origin_bam (str) :  path to extract region from to re-align

    target_bam(str) : path to write re-aligned reads to

    contig (str) : contig of selected region

    region_start (int) : start coordinate of region to re align (1 based)

    region_end (int) :end coordiante of  selected region (1 based)

    known_variants_vcf_path (str) : path to vcf containing reference variants

    interval_path (str) : Use this intervals to perform realignment, when not specified intervals are generated using RealignerTargetCreator

    interval_write_path (str) : when interval_path is not supplied, write the interval file here

    java_cmd (str) : Command to open java

    gatk_path (str) : path to GenomeAnalysisTK.jar

    reference (str) : path to reference Fasta file
    """
    if not os.path.exists('./gatk_tmp'):
        try:
            os.path.makedirs('./gatk_tmp')
        except Exception as e:
            pass

    # Detect reference from bam file
    if reference is None:
        reference = get_reference_path_from_bam(origin_bam)
    if reference is None:
        raise ValueError('Supply a path to a reference Fasta file (reference)')

    if interval_path is None:
        if interval_write_path is None:
            interval_write_path = f"{target_bam.replace('.bam','')}.intervals"

        target_creator_cmd = f'{java_cmd} {gatk_path} \
        -T RealignerTargetCreator \
        -R {reference} \
        -L {contig}:{region_start}-{region_end} \
        -known {known_variants_vcf_path} \
        -I {origin_bam} \
        -o {interval_write_path}'

        # Create the intervals file
        os.system(target_creator_cmd)
        interval_path = interval_write_path

    # Perform realignment
    realign_cmd = f'{java_cmd} {gatk_path} \
    -T IndelRealigner \
    -R {reference} \
    -targetIntervals {interval_path} \
    -known {known_variants_vcf_path} \
    -L {contig}:{region_start}-{region_end} \
    -I {origin_bam} \
    -o {target_bam} \
    -dcov 1000000 \
    -maxReads 2000000'
    os.system(realign_cmd)
    return target_bam


def add_readgroups_to_header(
        origin_bam_path,
        readgroups_in,
        target_bam_path=None,
        header_write_mode='auto'):
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
    if isinstance(readgroups_in, set):
        readGroupsDict = {}
        for readGroup in readgroups_in:
            flowCell, lane, sampleLib = readGroup.split('.')
            try:
                library, _ = sampleLib.rsplit('_', 1)
            except Exception as e:
                # the library is not part of the sample name:
                library = 'undefinedLibrary'
            readGroupsDict[readGroup] = {
                'ID': readGroup,
                'LB': library,
                'PL': 'ILLUMINA',
                'SM': sampleLib,
                'PU': readGroup}
    elif isinstance(readgroups_in, dict):
        readGroupsDict = readgroups_in
    else:
        raise ValueError("supply a set or dict for readgroups_in")

    # Write the re-headered bam file to this path
    complete_temp_path = origin_bam_path.replace('.bam', '') + '.rehead.bam'

    # When header_write_mode is auto, when samtools is available, samtools
    # will be used, otherwise pysam
    if header_write_mode == 'auto':
        if which('samtools') is None:
            header_write_mode = 'pysam'
        else:
            header_write_mode = 'samtools'

    if header_write_mode == 'pysam':
        with pysam.AlignmentFile(origin_bam_path, "rb") as origin:
            header = origin.header.copy()
            hCopy = header.to_dict()
            hCopy['RG'] = list(readGroupsDict.values())
            with pysam.AlignmentFile(complete_temp_path, "wb", header=hCopy) as out:
                for read in origin:
                    out.write(read)

        os.rename(complete_temp_path, target_bam_path)

    elif header_write_mode == 'samtools':
        with pysam.AlignmentFile(origin_bam_path, "rb") as origin:
            header = origin.header.copy()

            # Write the new header to this sam file:
            headerSamFilePath = origin_bam_path.replace(
                '.bam', '') + '.header.sam'

            hCopy = header.to_dict()
            hCopy['RG'] = list(readGroupsDict.values())

            # Write the sam file with the complete header:
            with pysam.AlignmentFile(headerSamFilePath, 'w', header=hCopy):
                pass

        # Concatenate and remove origin
        rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {origin_bam_path}; }} | samtools view -b > {complete_temp_path} ;
                mv {complete_temp_path } {target_bam_path};rm {headerSamFilePath};
                """
        os.system(rehead_cmd)
    else:
        raise ValueError(
            'header_write_mode should be either, auto, pysam or samtools')
