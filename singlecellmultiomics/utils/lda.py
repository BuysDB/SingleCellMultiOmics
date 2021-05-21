import pandas as pd
import numpy as np
from sklearn.decomposition import LatentDirichletAllocation

class SCMO_LDA(LatentDirichletAllocation):
    """
    This class is a slight expansion of the sklearn LatentDirichletAllocation,
    and implements two scmo specific additions:

    - reconstruction of the count matrix using reconstruct()
    - keep count matrix structure when performing fit_transform

    Example:
        >>> from singlecellmultiomics.utils import SCMO_LDA
        >>> SCMO_LDA(n_jobs=-1)
        >>> topics = lda.fit_transform(X)
        >>> imputed_X = lda.reconstruct(X)
    """

    def fit_transform(self, X: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """
        Expects a sample by gene/genomic location pd.dataframe
        """
        return pd.DataFrame(
                LatentDirichletAllocation.fit_transform(
                    self,
                    X,
                    **kwargs),
                index=X.index
                )

    def reconstruct(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Reconstruct imputed count matrix
        This method converts the LDA components to probabilities
        and then multiplies these by the total counts observed for
        each cell

        Example:
            >>> from singlecellmultiomics.utils import SCMO_LDA
            >>> SCMO_LDA(n_jobs=-1)
            >>> topics = lda.fit_transform(X)
            >>> imputed_X = lda.reconstruct(X)

        Warning:
            Make sure to call .fit_transform or .fit first!

        Args:
            X: sample by gene/genomic location pandas DataFrame

        Returns:
            reconstructed_count_frame : sample by gene/genomic location pandas DataFrame

        """
        # Obtain topic weights for each cell:
        tf = self.transform(X)
        # Convert lda components to probabilities
        comps = pd.DataFrame(
            self.components_ / self.components_.sum(axis=1)[:, np.newaxis],
            columns=X.columns )

        # Convert the probabilities to counts and reconstruct the count matrix
        return pd.DataFrame([
                (comps.T*row).sum(1)*total_cuts
                for (cell, row), total_cuts in
                    zip(pd.DataFrame(tf,index=X.index).iterrows(),
                    X.sum(1)
                )
            ], index=X.index )
