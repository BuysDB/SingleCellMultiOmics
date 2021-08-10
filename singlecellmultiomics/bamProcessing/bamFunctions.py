import os
import pysam
import time
import contextlib
from shutil import which, move
from singlecellmultiomics.utils import BlockZip, Prefetcher, get_file_type
import uuid
import os
from collections import defaultdict, Counter
from singlecellmultiomics.bamProcessing.pileup import pileup_truncated
import numpy as np
import pandas as pd
from typing import Generator
from multiprocessing import Pool


def get_index_path(bam_path: str):
    """
    Obtain path to bam index

    Returns:
        path_to_index(str) : path to the index file, None if not available
    """
    for p in [bam_path+'.bai', bam_path.replace('.bam','.bai')]:
        if os.path.exists(p):
            return p
    return None


def _get_r1_counts_per_cell(args):
    """Obtain the amount of unique read1 reads per cell (Function is used as Pool chunk)

    Args:
        args: bam_path, contig, prefix

    Returns:
        cell_obs (Counter) : {sampleA:n_molecules, sampleB:m_molecules, ...}
    """
    bam_path, contig, prefix = args
    cell_obs = Counter()


    with pysam.AlignmentFile(bam_path) as alignments:
        for read in alignments.fetch(contig):
            if read.is_qcfail or read.is_duplicate or not read.is_read1:
                continue
            if prefix is not None:
                cell_obs[prefix, read.get_tag('SM')]+=1
            else:
                cell_obs[read.get_tag('SM')]+=1
    return cell_obs

def get_r1_counts_per_cell(bam_path, prefix_with_bam=False):
    """Obtain the amount of unique read1 reads per cell

    Args:
        bam_path : str
        prefix_with_bam(bool) : add bam name as prefix of cell name
    Returns:
        cell_obs (Counter) : {sampleA:n_molecules, sampleB:m_molecules, ...}
    """

    if type(bam_path)==str:
        bam_paths = [bam_path]
    else:
        bam_paths=bam_path


    cell_obs = Counter()

    def generate_commands(bam_paths, prefix_with_bam):

        for bam_path in bam_paths:
            if prefix_with_bam:
                prefix = bam_path.split('/')[-1].replace('.bam','')
            else:
                prefix=None

            for contig in get_contigs_with_reads(bam_path):
                yield (bam_path, contig, prefix)



    with Pool() as workers:
        for cell_obs_for_contig in workers.imap_unordered(
            _get_r1_counts_per_cell,
            generate_commands(bam_paths, prefix_with_bam)):

            cell_obs += cell_obs_for_contig

    return cell_obs


def get_contigs(bam: str)  -> Generator :
    """
    Get all contigs listed in the bam file header

    Args:
        bam_path(str): path to bam file or pysam handle


    Returns:
        contigs(str) list

    """

    if type(bam) is str:
        with pysam.AlignmentFile(bam) as a:
            references = a.references
        return references
    elif type(bam) is pysam.AlignmentFile:
        return bam.references

    raise ValueError('Supply either path to bam or pysam.AlignmentFile')




def get_contigs_with_reads(bam_path: str, with_length: bool = False)  -> Generator :
    """
    Get all contigs with reads mapped to them

    Args:
        bam_path(str): path to bam file

        with_length(bool): also yield the length of the contig

    Yields:
        contig(str)

    """
    for line in pysam.idxstats(bam_path).split('\n'):
        try:
            contig, contig_len, mapped_reads, unmapped_reads = line.strip().split()
            mapped_reads, unmapped_reads = int(mapped_reads), int(unmapped_reads)
            if mapped_reads>0 or unmapped_reads>0:
                if with_length:
                    yield contig, int(contig_len)
                else:
                    yield contig
        except ValueError:
            pass

def merge_bams( bams: list, output_path: str, threads: int=4 ):
    """Merge bamfiles to output_path

    When a single bam file is supplied, the bam file is moved to  output_path
    All input bam files are removed

    Args:
        bams : list or tuple containing paths to bam files to merge
        output_path (str): target path

    Returns:
        output_path (str)

    """
    assert threads>=1
    if len(bams) == 1:
        assert os.path.exists(bams[0]+'.bai'), 'Only indexed files can be merged'
        move(bams[0], output_path)
        move(bams[0]+'.bai', output_path+'.bai')
    else:
        assert all((os.path.exists(bams[0]+'.bai') for bam in bams)), 'Only indexed files can be merged'
        if which('samtools') is None:
            pysam.merge(output_path, *bams, f'-@ {threads} -f -p -c') #-c to only keep the same id once
        else:
            # This above command can have issues...
            os.system(f'samtools merge {output_path} {" ".join(bams)} -@ {threads} -f -p -c')

        pysam.index(output_path, f'-@ {threads}')
        for o in bams:
            os.remove(o)
            os.remove(o+'.bai')
    return output_path

def verify_and_fix_bam(bam_path):
    """
    Check if the bam file is not truncated and indexed.
    Also regenerates index when its older than the bam file
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

        # Check index modification time
        if not index_missing:
            index_path = get_index_path(bam_path)
            if index_path is None or os.path.getmtime(index_path) < os.path.getmtime(bam_path):
                index_missing = True

        if index_missing:
            print(f"The bam file {bam_path} has no or an older index, attempting an index re-build ..")
            pysam.index(bam_path)
            print(f"Succes!")

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
    try:
        return set([entry['SM'] for entry in handle.header.as_dict()['RG']])
    except Exception as e:
        samples = set()
        for read in handle:
            samples.add( read.get_tag('SM') )
        return samples

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

def get_read_group_to_sample_dict(bam):
    """ Obtain a dictionary containing {'read_group' : 'sample' , ...}
        Args:
            bam_file(pysam.AlignmentFile) or path to bam file or pysam object
    """
    d = get_sample_to_read_group_dict(bam)
    r2s = {}
    for sample, read_groups in d.items():
        for rg in read_groups:
            r2s[rg] = sample

    return r2s

def get_read_group_format(bam):
    """Obtain read group format

    Args:
        bam (str) : path to bam file

    Returns:
        type (int) : read group format

    Raises:
        ValueError : when read group format cannot be determined
    """
    with pysam.AlignmentFile(bam) as f:
        d = f.header.as_dict()
        format_type_obs = Counter()
        for read_group_dict in d['RG']:

            if read_group_dict.get('LB','')==read_group_dict.get('SM',''):
                format_type_obs[1]+=1
            else:
                format_type_obs[0]+=1

        if len(format_type_obs)==1:
            return format_type_obs.most_common(1)[0][0]
        raise ValueError('Failed to indentify read group format ')


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


def get_contig_size(bam, contig):
    """Extract the length of a contig from a bam file

    Args:
        bam (str or pysam.AlignmentFile) : handle to bam file or path to bam file
        contig (str)

    Returns:
        length (int)
    """
    if type(bam) is str:
        with pysam.AlignmentFile(bam) as a:
            size = get_contig_size(a, contig)
        return size
    elif type(bam) is pysam.AlignmentFile:
        for c,length in zip(bam.references, bam.lengths):
            if c==contig:
                return length
    else:
        raise ValueError(
            'Supply either a path to a bam file or pysam.AlignmentFile object')

    return None

def get_contig_sizes(bam):
    """Extract lengths of all contigs from a bam file

    Args:
        bam (str or pysam.AlignmentFile) : handle to bam file or path to bam file

    Returns:
        contig_lengths : dict (contig:length (int) )
    """

    if type(bam) is str:
        with pysam.AlignmentFile(bam) as a:
            sizes = get_contig_sizes(a)
        return sizes
    elif type(bam) is pysam.AlignmentFile:
        return dict(zip(bam.references, bam.lengths))
    else:
        raise ValueError(
            'Supply either a path to a bam file or pysam.AlignmentFile object')


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
        local_temp_sort=True,
        input_is_sorted=False,
        mode='wb',
        fast_compression=False, # Use fast compression for merge (-1 flag)
        temp_prefix = 'SCMO',
        **kwargs
        ):
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

        input_is_sorted(bool) : Assume the input is sorted, no sorting will be applied

        mode (str) : Output mode, use wbu for uncompressed writing.

        **kwargs : arguments to pass to the new pysam.AlignmentFile output handle

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
    elif type(origin_bam) is str:
        with pysam.AlignmentFile(origin_bam) as ain:
            header = ain.header.copy()
    elif origin_bam is not None:
        header = origin_bam.header.copy()
    else:
        raise ValueError("Supply a header or origin_bam object")


    try: # Try to open the temp unsorted file:
        unsorted_alignments = pysam.AlignmentFile(
            unsorted_path, mode, header=header, **kwargs)
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
    if input_is_sorted is False:
        sort_and_index(
            unsorted_path,
            write_path,
            remove_unsorted=True,
            local_temp_sort=local_temp_sort,
            fast_compression=fast_compression,
            prefix=temp_prefix
            )
    else:
        os.rename(unsorted_path, write_path)

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

    blocked = set()
    for prog_entry in input_header['PG']:
        if prog_entry.get('PN','') == program_name:
            blocked.add( prog_entry.get('ID','') )

    proposed_id = program_name
    i = 0
    while proposed_id in blocked:
        proposed_id= f'{program_name}_{i}'
        i+=1

    input_header['PG'].append({
        'ID': proposed_id,
        'PN': program_name,
        'CL': command_line,
        'VN': version,
        'DS': description
    })

def  add_blacklisted_region(input_header, contig=None, start=None, end=None, source='auto_blacklisted'):
    if 'CO' not in input_header: # Blacklisted regions, as freetext comment
        input_header['CO'] = []
    input_header['CO'].append(f'scmo_blacklisted\t{contig}:{start} {end}\tsource:{source}')

def get_blacklisted_regions_from_bam(bam_path):
    blacklisted = {}
    with pysam.AlignmentFile(bam_path) as alignments:
        if 'CO' in alignments.header:
            for record in alignments.header['CO']:
                parts = record.strip().split('\t')
                if parts[0]!='scmo_blacklisted':
                    continue

                region = parts[1]
                contig, span = region.split(' ')
                start, end = span.split('-')
                start = int(start)
                end = int(end)
                if not contig in blacklisted:
                    blacklisted[contig] = []
                blacklisted[contig].append( (start,end) )
    return blacklisted

def bam_is_processed_by_program(alignments, program='bamtagmultiome'):
    """Check if bam file has been processed by the supplied program

    This function checks if there is an entry available in the 'PG'
    block in header of the bam file where PN matches the program supplied.

    Args:
        alignments(pysam.AlignmentFile) : handle to bam file
        program(str) : program to look for

    Returns:
        program_present(bool) : the program was used to process the bam file

    """
    for program_dict in alignments.header.as_dict()['PG']:
        if program_dict.get('PN','') ==  program:
            return True
    return False


def sort_and_index(
        unsorted_path,
        sorted_path,
        remove_unsorted=False,
        local_temp_sort=True,
        fast_compression=False,
        prefix='TMP'
        ):
    """ Sort and index a bam file
    Args:
        unsorted_path (str) : path to unsorted bam file

        sorted_path (str) : write sorted file here

        remove_unsorted (bool) : remove the unsorted file

        local_temp_sort(bool): create temporary files in target directory
    Raises:
        SamtoolsError when sorting or indexing fails
    """

    if prefix is None:
        raise 'prefix cannot be None'

    if local_temp_sort:
        base_directory = os.path.abspath( os.path.dirname(sorted_path) )
        prefix = f'{prefix}.{uuid.uuid4()}'
        temp_path_first = f'{base_directory}/{prefix}'
        if temp_path_first.startswith('/TMP'):
            # Perform sort in current directory
            temp_path_first = f'./{prefix}'


        # Try to sort at multiple locations, if sorting fails try the next until all locations have been tried
        temp_paths =  [temp_path_first, f'/tmp/{prefix}', f'./{prefix}' ]
        for i, temp_path in enumerate(temp_paths):
            failed = False
            try:
                pysam.sort(
                    '-o',
                    sorted_path,
                    '-T',
                    f'{temp_path}', # Current directory with a random prefix
                    unsorted_path, ('-l 1' if fast_compression else '-l 3')
                )
            except Exception as e:
                print(f'Sorting failed at temp: {temp_path}')
                failed = True
                if i==len(temp_paths)-1:
                    raise

            if not failed:
                break
    else:
        pysam.sort("-o", sorted_path, unsorted_path, ('-l 1' if fast_compression else '-l 3'))
    pysam.index(sorted_path, '-@ 4')
    if remove_unsorted:
        os.remove(unsorted_path)


class MapabilityReader(Prefetcher):

    def __init__(self, mapability_safe_file_path, read_all=False, dont_open=True):
        self.args = locals().copy()
        self.handle = None
        del self.args['self']
        self.mapability_safe_file_path = mapability_safe_file_path
        if not dont_open:
            self.handle = BlockZip(mapability_safe_file_path, 'r')


    def instance(self, arg_update=None):

        if 'self' in self.args:
            del self.args['self']

        if arg_update is not None:
            self.args.update(arg_update)

        clone = MapabilityReader(**self.args, dont_open=False)
        return clone

        # Todo: exit statements
    def prefetch(self, contig, start, end):
        clone = self.instance()
        clone.handle.read_contig_to_cache( contig, region_start=start, region_end=end)
        return clone


    def __getitem__(self, contig_ds_strand):
        if self.handle is None:
            self.handle = BlockZip(self.mapability_safe_file_path, 'r')

        contig, ds, strand = contig_ds_strand
        return self.handle[contig, ds, strand]

    def site_is_mapable(self, contig, ds, strand):

        if self.handle is None:
            self.handle = BlockZip(self.mapability_safe_file_path, 'r')

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

def get_random_locations(bam, n):
    """Select random locations in the supplied bam file

    bam(str or pysam.AlignmentFile)

    n(int) : amount of locations to generate

    returns: generator of (contig,position) tuples
    """

    cs = get_contig_sizes(bam)
    # Obtain cumulative amount of bases per contig:
    cumulative_size = np.cumsum([size for contig,size in cs.items()])
    cs_contigs = np.array( list(cs.keys()) )

    random_locations = np.random.randint(0, cumulative_size[-1], n)
    indices = np.searchsorted(cumulative_size, random_locations)
    random_positions = random_locations-np.concatenate( ([0], cumulative_size))[indices]
    random_contigs = cs_contigs[indices]

    return zip(random_contigs, random_positions)

def sam_to_bam(sam_in, bam_out, threads = 4):
    """
    Convert sam file to sorted bam file

    Args:

        sam_in(str) : input sam file path

        bam_out(str) : output bam file path

    """
    os.system(f'samtools view {sam_in} -b | samtools sort -@ {threads} > {bam_out}; samtools index {bam_out};')


def sample_location(handle, contig, pos,dedup=True, qc=True):
    """
    Obtain dictionary containing the coverage for every sample

    Args:
        handle (pysam.AlignmentFile)  : File to obtain reads from

        contig (str) : contig to sample

        pos (int) : coordinate to sample (zero based)

        dedup(bool) : ignore duplicated reads

        qc(bool) : ignore qc failed reads

    Returns:
        samples (dict) : dictionary with amount of reads at the selected locations

    """
    overlap = {} # sample -> coverage
    for column in pileup_truncated(handle, contig, pos, pos+1):
        for pileread in column.pileups:
            if not pileread.is_del and not pileread.is_refskip:
                sample = pileread.alignment.get_tag('SM')
                if dedup and pileread.alignment.is_duplicate:
                    continue
                if qc and pileread.alignment.is_qcfail:
                    continue

                if not sample in overlap:
                    overlap[sample]=1
                else:
                    overlap[sample]+=1
    return overlap


def get_read_group_from_read(read, format, with_attr_dict=False):
    rg_id=None
    platform = read.get_tag('Fc') if read.has_tag('Fc') else 'NONE'
    lane = read.get_tag('La') if read.has_tag('La') else 'NONE'
    library = read.get_tag('LY') if read.has_tag('LY') else 'NONE'
    if format==0:
        sample = read.get_tag('SM')
    elif format==1:
        sample = read.get_tag('LY')
    else:
        raise ValueError(f'{format} is an unknown read group format')

    rg_id =  f'{platform}.{lane}.{sample}'
    if with_attr_dict:
        rg_dict = {
                'ID':rg_id,
                'LB':library,
                'PL':platform,
                'SM':sample,
                'PU': rg_id }

    if rg_id is None:
        raise ValueError('No read group information could be extracted')

    if with_attr_dict:
        return rg_id, rg_dict
    return rg_id


def random_sample_bam(bam,n,**sample_location_args):
    """Sample a bam file at random locations

    Args:
        bam (pysam.AlignmentFile) : bam to sample from

        n (int) : Amount of locations to sample

        *sample_location_args : arguments to pass to sample_location

    Returns:
        samples (dict) : dictionary with amount of reads at the selected locations

    """
    r = [sample_location(bam, contig, pos,**sample_location_args)
    for contig, pos in get_random_locations(bam, n)]
    return pd.DataFrame(r)


def replace_bam_header(origin_bam_path, header, target_bam_path=None, header_write_mode='auto'):

    if target_bam_path is None:
        target_bam_path = origin_bam_path

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

        with pysam.AlignmentFile(complete_temp_path, "wb", header=header) as out, \
             pysam.AlignmentFile(origin_bam_path) as origin:
            for read in origin:
                out.write(read)

        os.rename(complete_temp_path, target_bam_path)

    elif header_write_mode == 'samtools':

        # Write the new header to this sam file:
        headerSamFilePath = origin_bam_path.replace(
            '.bam', '') + '.header.sam'

        # Write the sam file with the complete header:
        with pysam.AlignmentFile(headerSamFilePath, 'w', header=header):
            pass

        # Concatenate and remove origin
        rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {origin_bam_path}; }} | samtools view -b > {complete_temp_path} ;
                mv {complete_temp_path } {target_bam_path};rm {headerSamFilePath};
                """
        os.system(rehead_cmd)
    else:
        raise ValueError(
            'header_write_mode should be either, auto, pysam or samtools')



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

    with pysam.AlignmentFile(origin_bam_path, "rb") as origin:
        header = origin.header.copy()
        hCopy = header.to_dict()
        hCopy['RG'] = list(readGroupsDict.values())
        replace_bam_header(origin_bam_path, hCopy,
                            target_bam_path=target_bam_path,
                            header_write_mode=header_write_mode)



def mate_iter(alignments, **kwargs):
    buffer =  dict()
    for read in alignments.fetch(**kwargs):

        if not read.is_paired or not read.is_proper_pair:
            yield [read, None]

        if read.is_qcfail or read.is_supplementary:
            continue
            #yield read

        if read.query_name in buffer:
            if read.is_read1:
                yield read, buffer.pop(read.query_name)
            elif read.is_read2:
                yield buffer.pop(read.query_name), read
        else:
            buffer[read.query_name] = read
