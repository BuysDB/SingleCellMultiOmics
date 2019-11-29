def molecule_iterator_filter(
    molecule_iterator,
    min_fragments=None,
    max_fragments=None,
    min_mapping_qual=None,
    max_mapping_qual=None,
    min_ivt_duplicates=None,
    max_ivt_duplicates=None,
    both_pairs_mapped=None,
    min_span=None,
    max_span=None
):
    """ Filter iterable with molecules

    molecule_iterator (iterable) : molecules to filter from
    min_fragments (int) : minimum amount of fragments associated to molecule
    max_fragments (int) : maximum amount of fragments associated to molecule
    min_mapping_qual(int) : minimum maximum mapping quality for a single associated fragment
    max_mapping_qual(int) : maximum maximum mapping quality for a single associated fragment
    min_ivt_duplicates(int) : minimum amount of in vitro transcription copies
    max_ivt_duplicates(int) : maximum amount of in vitro transcription copies
    both_pairs_mapped(bool) : molecule should have at least one fragment with both pairs mapped
    min_span(int) : minimum amount of bases aligned with reference
    max_span : maximum amount of bases aligned with reference

    """

    for molecule in molecule_iterator:
        if min_fragments is not None and len(molecule)<min_fragments:
            continue
        if max_fragments is not None and len(molecule)>max_fragments:
            continue

        if min_mapping_qual is not None and molecule.get_max_mapping_qual()<min_mapping_qual:
            continue
        if max_mapping_qual is not None and molecule.get_max_mapping_qual()>max_mapping_qual:
            continue


        if min_ivt_duplicates is not None and len( molecule.get_rt_reactions() ) <  min_ivt_duplicates:
            continue
        if min_ivt_duplicates is not None and len( molecule.get_rt_reactions() ) >  max_ivt_duplicates:
            continue

        if both_pairs_mapped is not None:
            found_both = False
            for fragment in molecule:
                if fragment.has_R1() and fragment.has_R2():
                    found_both=True
            if not found_both:
                continue

        if min_span is not None and molecule.get_safely_aligned_length()<min_span:
            continue

        if max_span is not None and molecule.get_safely_aligned_length()>max_span:
            continue

        yield molecule
