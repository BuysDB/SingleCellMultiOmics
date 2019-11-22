import sklearn.ensemble
import numpy as np


def get_consensus_training_data(
        molecule_iterator,
        mask_variants=None,
        n_train=100_000,
        **feature_matrix_args):
    X = None
    y = []

    last_end = None
    last_chrom = None
    for i, molecule in enumerate(molecule_iterator):
        # Never train the same genomic location twice
        if last_chrom is not None and last_chrom != molecule.chromosome:
            last_end = None
        if last_end is not None and molecule.spanStart < last_end:
            continue
        x, _y = molecule.get_base_calling_training_data(
            mask_variants, **feature_matrix_args)
        if X is None:
            X = np.empty((0, x.shape[1]))
            print(
                f"Creating feature matrix with {x.shape[1]} dimensions and {n_train} training base-calls")
        y += _y
        X = np.append(X, x, axis=0)
        last_chrom = molecule.chromosome
        if molecule.spanEnd is not None:
            last_end = molecule.spanEnd
        else:
            last_end += len(_y)
        if len(X) >= n_train:
            break
    print(
        f'Finished, last genomic coordinate: {molecule.chromosome} {molecule.spanEnd}')
    return X, y


def train_consensus_model(
        molecule_iterator,
        mask_variants=None,
        classifier=None,
        n_train=100_000,
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
        molecule_iterator, mask_variants=mask_variants, n_train=n_train, **feature_matrix_args)
    y = np.array(y)
    # remove unkown ref bases from set
    X = np.array(X)[y != 'N']
    y = y[y != 'N']
    classifier.fit(X, y)
    if isinstance(classifier, sklearn.ensemble.forest.RandomForestClassifier):
        print(f"Model out of bag accuracy: {classifier.oob_score_}")
    classifier.n_jobs = 1  # fix amount of jobs to one, otherwise apply will be very slow
    return classifier
