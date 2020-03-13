from singlecellmultiomics.bamProcessing.pileup import pileup_truncated

def get_pileup_vect(alignments, contig, pos, ref, alt):
    """Create feature vector for selected variant

     Args:
        alignments(pysam.AlignmentFile) : Handle to alignmentfile
        contig(str) : contig to perform pileup
        pos(int) : zeros based position of variant to pileup
        ref(str) : reference base
        alt(str) : alternative base

    Returns
        total(int) : Total amount of bases overlapping with selected location
        ref_calls  : Total amount of bases matching ref
        alt_calls : Total amount of bases matching alt
        other_calls : Total amount of bases matching neither ref nor alt
    """
    total = 0
    ref_calls = 0
    alt_calls = 0
    other_calls = 0
    start=pos
    stop = pos+1
    for pileupcolumn in pileup_truncated(alignments,contig,start,stop,stepper='all'):
        for i,pileupread in enumerate(pileupcolumn.pileups):
            if not pileupread.is_del and not pileupread.is_refskip:
                call = pileupread.alignment.query_sequence[pileupread.query_position]
                if call==ref:
                    ref_calls += 1
                elif call==alt:
                    alt_calls += 1
                else:
                    other_calls += 1
                other_calls+=1
    return total, ref_calls, alt_calls, other_calls


# Create mapping quality feature vector:
def get_mapping_q_vect(alignments_handle, contig, pos, radius=150):

    """Obtain histogram of mapping qualties, clipped at 60

     Args:
        alignments(pysam.AlignmentFile) : Handle to alignmentfile
        contig(str) : contig
        pos(int) : zeros based position of location to check mapping qualties
        radius(int) : radius to check around selected location

    Returns:
        mapping_qualities(list) : Histogram with 7 bins (0 to highest mapping quality)
    """
    mapping_qualities = [0]*7
    for read in alignments_handle.fetch(contig, pos-radius, pos+radius):
        mapping_qualities[min(60,int(read.mapping_quality/10))]+=1
    return mapping_qualities
