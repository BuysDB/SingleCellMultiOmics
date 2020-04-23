import sklearn.ensemble
import numpy as np
import pandas as pd
import time

def calculate_consensus(molecule, consensus_model, molecular_identifier, out, **model_kwargs ):
    """
    Create consensus read for molecule

    Args:
        molecule (singlecellmultiomics.molecule.Molecule)

        consensus_model

        molecular_identifier (str) : identier for this molecule, will be suffixed to the reference_id

        out(pysam.AlingmentFile) : target bam file

        **model_kwargs : arguments passed to the consensus model

    """
    try:
        consensus_reads = molecule.deduplicate_to_single_CIGAR_spaced(
            out,
            f'c_{molecule.get_a_reference_id()}_{molecular_identifier}',
            consensus_model,
            NUC_RADIUS=model_kwargs['consensus_k_rad']
            )
        for consensus_read in consensus_reads:
            consensus_read.set_tag('RG', molecule[0].get_read_group())
            consensus_read.set_tag('mi', molecular_identifier)
            out.write(consensus_read)
    except Exception as e:

        molecule.set_rejection_reason('CONSENSUS_FAILED',set_qcfail=True)
        molecule.write_pysam(out)

def base_calling_matrix_to_df(
        x,
        ref_info=None,
        NUC_RADIUS=1,
        USE_RT=True):
    """
    Convert numpy base calling feature matrix to pandas dataframe with annotated columns

    Args:
        x(np.array) : feature matrix
        ref_info(list) : reference position annotations (will be used as index)
        NUC_RADIUS(int) : generate kmer features target nucleotide
        USE_RT(bool) : use RT reaction features

    Returns:
        df (pd.DataFrame)
    """
    df = pd.DataFrame(x)
    # annotate the columns
    BASE_COUNT = 5
    RT_INDEX = 7 if USE_RT else None
    STRAND_INDEX = 0
    PHRED_INDEX = 1
    RC_INDEX = 2
    MATE_INDEX = 3
    CYCLE_INDEX = 4
    MQ_INDEX = 5
    FS_INDEX = 6
    COLUMN_OFFSET = 0
    features_per_block = 8 - (not USE_RT)

    block_header = ["?"] * features_per_block
    block_header[STRAND_INDEX] = 'strand'
    block_header[PHRED_INDEX] = 'phred'
    block_header[RC_INDEX] = 'read_count'
    block_header[MATE_INDEX] = 'mate'
    block_header[CYCLE_INDEX] = 'cycle'
    block_header[MQ_INDEX] = 'mq'
    block_header[FS_INDEX] = 'fragment_size'
    block_header[RT_INDEX] = 'rt_reactions'
    k_header = []
    for k in range(NUC_RADIUS * 2 + 1):
        for base in 'ACGTN':
            k_header += [(k, b, base) for b in block_header]

    try:
        df.columns = pd.MultiIndex.from_tuples(k_header)
    except ValueError:  # the dataframe is a concateenation of multiple molecules
        pass
    if ref_info is not None:
        df.index = pd.MultiIndex.from_tuples(ref_info)
    return df

def get_consensus_training_data(
        molecule_iterator,
        mask_variants=None,
        n_train=100_000, # When None, training data is created until molecule source depletion
        skip_already_covered_bases = True,
        #yield_results=False, # Yield results instead of returning a matrix
        **feature_matrix_args):

    """
    Create a tensor/matrix containing alignment and base calling information, which can be used for consensus calling.
    This function also creates a vector containing the corresponding reference bases, which can be used for training a consensus model.

    Args:
        molecule_iterator : generator which generates molecules from which base calling feature matrices are extracted

        mask_variants (pysam.VariantFile) : variant locations which should be excluded from the matrix

        n_train (int) : amount of rows in the matrix

        skip_already_covered_bases(bool) : when True every reference position is at most a single row in the output matrix, this prevents overfitting

        **feature_matrix_args : Arguments to pass to the feature matrix function of the molecules.

    """

    #if not yield_results:
    X = None
    y = []
    molecules_used = 0
    training_set_size = 0
    last_end = None
    last_chrom = None
    try:

        for i, molecule in enumerate(molecule_iterator):
            # Never train the same genomic location twice
            if skip_already_covered_bases:
                if last_chrom is not None and last_chrom != molecule.chromosome:
                    last_end = None
                if last_end is not None and molecule.spanStart < last_end:
                    continue

            train_result = molecule.get_base_calling_training_data(
                mask_variants, **feature_matrix_args)

            if train_result is None:
                # Continue when the molecule does not have bases where we can learn from
                continue

            x, _y = train_result
            training_set_size += len(_y)

            #if yield_results:
            #    yield x, _y

            #else:
            if X is None:
                X = np.empty((0, x.shape[1]))
                print(
                    f"Creating feature matrix with {x.shape[1]} dimensions and {n_train} training base-calls")
            y += _y
            X = np.append(X, x, axis=0)
            last_chrom = molecule.chromosome

            if training_set_size >= n_train:
                break

            molecules_used+=1
            if molecule.spanEnd is not None:
                last_end = molecule.spanEnd
            else:
                last_end += len(_y)


    except KeyboardInterrupt as e:
        print("Got keyboard interrupt, stopping to load more data")

    if molecules_used > 0:
        print(
                f'Finished, last genomic coordinate: {molecule.chromosome} {molecule.spanEnd}, training set size is {training_set_size}, used {molecules_used} molecules for training')
    #if not yield_results:
    return X, y





def train_consensus_model(
        molecule_iterator,
        mask_variants=None,
        classifier=None,
        n_train=100_000,
        skip_already_covered_bases = True,
        **feature_matrix_args):
    if classifier is None:  # default to random forest
        classifier = sklearn.ensemble.RandomForestClassifier(
            n_jobs=-1,
            n_estimators=100,
            oob_score=True,
            max_depth=7,
            min_samples_leaf=5
        )
    X, y = get_consensus_training_data(
        molecule_iterator, mask_variants=mask_variants, n_train=n_train,
        skip_already_covered_bases=skip_already_covered_bases,**feature_matrix_args)
    y = np.array(y)
    # remove unkown ref bases from set
    X = np.array(X)[y != 'N']
    y = y[y != 'N']
    classifier.fit(X, y)
    if isinstance(classifier, sklearn.ensemble.forest.RandomForestClassifier):
        print(f"Model out of bag accuracy: {classifier.oob_score_}")
    classifier.n_jobs = 1  # fix amount of jobs to one, otherwise apply will be very slow
    return classifier
