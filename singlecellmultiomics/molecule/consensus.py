import sklearn.ensemble
import numpy as np

def get_consensus_training_data(molecule_iterator, mask_variants=None, n_train=100_000):
    X = None
    y = []

    for i,molecule in enumerate( molecule_iterator ):
        x,_y = molecule.get_base_calling_training_data(mask_variants)
        if X is None:
            X = np.empty((0,x.shape[1]))
            print(f"Creating feature matrix with {x.shape[1]} dimensions and {n_train} training base-calls")
        y+=_y
        X = np.append(X,x,axis=0)
        if len(X)>=n_train:
            break
    return X,y

def train_consensus_model(molecule_iterator, mask_variants=None, classifier=None, n_train=100_000):
    if classifier is None: # default to random forest
        classifier = sklearn.ensemble.RandomForestClassifier(
                n_jobs=8,
                n_estimators=100,
                oob_score=True
                )
    X,y = get_consensus_training_data(molecule_iterator, mask_variants=mask_variants, n_train=n_train)
    classifier.fit(X,y)
    if type(classifier)==sklearn.ensemble.forest.RandomForestClassifier:
        print(f"Model out of bag accuracy: {classifier.oob_score_}")
    return classifier
