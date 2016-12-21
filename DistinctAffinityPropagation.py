import numpy as np, pandas as pd, scipy
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from scipy.spatial.distance import euclidean


def affinity_propagation(X, S, max_cutoff, min_cutoff, preference=None, convergence_iter=15, max_iter=200,
                         damping=0.5, copy=True, verbose=False,
                         return_n_iter=False):
    """Perform Affinity Propagation Clustering of data

    Read more in the :ref:`User Guide <affinity_propagation>`.

    Parameters
    ----------

    S : array-like, shape (n_samples, n_samples)
        Matrix of similarities between points

    max_cutoff: the min cut_off to group sample together by affinity metrics

    min_cutoff: the max cut_off to find the most distinct pair by affinity metrics

    preference : array-like, shape (n_samples,) or float, optional
        Preferences for each point - points with larger values of
        preferences are more likely to be chosen as exemplars. The number of
        exemplars, i.e. of clusters, is influenced by the input preferences
        value. If the preferences are not passed as arguments, they will be
        set to the median of the input similarities (resulting in a moderate
        number of clusters). For a smaller amount of clusters, this can be set
        to the minimum value of the similarities.

    convergence_iter : int, optional, default: 15
        Number of iterations with no change in the number
        of estimated clusters that stops the convergence.

    max_iter : int, optional, default: 200
        Maximum number of iterations

    damping : float, optional, default: 0.5
        Damping factor between 0.5 and 1.

    copy : boolean, optional, default: True
        If copy is False, the affinity matrix is modified inplace by the
        algorithm, for memory efficiency

    verbose : boolean, optional, default: False
        The verbosity level

    return_n_iter : bool, default False
        Whether or not to return the number of iterations.

    Returns
    -------

    cluster_centers_indices : array, shape (n_clusters,)
        index of clusters centers

    labels : array, shape (n_samples,)
        cluster labels for each point

    n_iter : int
        number of iterations run. Returned only if `return_n_iter` is
        set to True.

    Notes
    -----
    This function should only be called by .fit function

    References
    ----------
    """

    n_samples = S.shape[0]
    cluster_centers_indices = []
    labels = []
    frontiers = np.zeros(n_samples)
    n_iter = 0

    affinity_matrix = S.copy()

    if S.shape[0] != S.shape[1]:
        raise ValueError("S must be a square array (shape=%s)" % repr(S.shape))
    # if preference is not None:
    #     first_distinct = preference
    # if damping < 0.5 or damping >= 1:
    #     raise ValueError('damping must be >= 0.5 and < 1')

    # find the most distinct pattern, if no distinct pattern is found, then return as one pattern
    while True and n_iter < max_iter:
        distinct_index = np.argmin(affinity_matrix)
        min_distance = affinity_matrix.flat[distinct_index]

        if min_distance > min_cutoff:
            break

        distinct_pairs = np.asarray(np.where(affinity_matrix == min_distance)).T
        distinct_set = remove_duplicate(distinct_pairs)

        distinct_set = filter_close_pair(S, cluster_centers_indices, distinct_set, min_cutoff)

        if len(distinct_set) == 0:
            break
        elif len(distinct_set) == 1:
            first_distinct = distinct_set[0]
        else:
            print "more than one pair are very different!"
            first_distinct = most_different_pair(X, distinct_set)

        cluster_centers_indices += first_distinct

        affinity_x = affinity_matrix[first_distinct[0], :]
        affinity_y = affinity_matrix[first_distinct[1], :]

        cluster_x = np.intersect1d(np.where(affinity_x > max_cutoff), np.where(affinity_x <=1))
        cluster_y = np.intersect1d(np.where(affinity_y > max_cutoff), np.where(affinity_y <=1))

        labels.append(cluster_x)
        labels.append(cluster_y)

        frontiers[cluster_x] = 1
        frontiers[cluster_y] = 1

        affinity_matrix[cluster_x, :] = float("inf")
        affinity_matrix[:, cluster_x] = float("inf")

        affinity_matrix[cluster_y, :] = float("inf")
        affinity_matrix[:, cluster_y] = float("inf")

        if np.sum(frontiers) == n_samples:
            break
        n_iter += 1
    return cluster_centers_indices, labels, n_iter

def remove_duplicate(pairs):
    """Remove duplicate pairs and return distinct dissimilar pairs

        Parameters
        ----------
        pairs: ndarray, contains pair of sample index.

        Return
        ----------
        a list of tuple containing pair of sample index.

        Notes
        -----
        This function should only be called by .fit function

        References
        ----------
        """

    distinct_set = set()

    for i in range(pairs.shape[0]):
        x, y = pairs[i]
        if (x, y) not in distinct_set and (y, x) not in distinct_set:
            distinct_set.add((x, y))
    distinct_set = list(distinct_set)
    return distinct_set

def most_different_pair(X, pairs):
    """Get the most distinct pair after the screening of first affinity metrics

            Parameters
            ----------
            X: ndarray, contains original sample data
            pairs: a list of pairs get from the first affinity metrics screening

            Return
            ----------
            a most distinct pair, [x, y]

            Notes
            -----

            References
            ----------
            """
    max_distance = float("-inf")
    distinct_pair = None
    for pair in pairs:
        sample_x, sample_y = X[pair[0], :], X[pair[1], :]
        cur_distance = euclidean(sample_x, sample_y)
        if cur_distance > max_distance:
            max_distance = cur_distance
            distinct_pair = pair
    return distinct_pair

def filter_close_pair(affinity_matrix, cluster_centers_indices, distinct_set, min_cutoff):
    """Make sure the new cluster are also very different with the old ones.

                Parameters
                ----------
                affinity_matrix: ndarray, the similarity matrix get from certain similarity function
                cluster_centers_indices: existing clusters index
                distinct_set: potential new cluster index
                min_cutoff: minimum similarity cutoff

                Return
                ----------
                a filtered distinct pair, [x, y]

                Notes
                -----

                References
                ----------
                """


###############################################################################

class DistinctAffinityPropagation():
    """Perform Affinity Propagation Clustering of data.

    Read more in the :ref:`User Guide <affinity_propagation>`.

    Parameters
    ----------
    damping : float, optional, default: 0.5
        Damping factor between 0.5 and 1.

    convergence_iter : int, optional, default: 15
        Number of iterations with no change in the number
        of estimated clusters that stops the convergence.

    max_iter : int, optional, default: 200
        Maximum number of iterations.

    copy : boolean, optional, default: True
        Make a copy of input data.

    preference : array-like, shape (n_samples,) or float, optional
        Preferences for each point - points with larger values of
        preferences are more likely to be chosen as exemplars. The number
        of exemplars, ie of clusters, is influenced by the input
        preferences value. If the preferences are not passed as arguments,
        they will be set to the median of the input similarities.

    affinity : string, optional, default=``euclidean``
        Which affinity to use. At the moment ``precomputed`` and
        ``euclidean`` are supported. ``euclidean`` uses the
        negative squared euclidean distance between points.

    verbose : boolean, optional, default: False
        Whether to be verbose.


    Attributes
    ----------
    cluster_centers_indices_ : array, shape (n_clusters,)
        Indices of cluster centers

    cluster_centers_ : array, shape (n_clusters, n_features)
        Cluster centers (if affinity != ``precomputed``).

    labels_ : array, shape (n_samples,)
        Labels of each point

    affinity_matrix_ : array, shape (n_samples, n_samples)
        Stores the affinity matrix used in ``fit``.

    n_iter_ : int
        Number of iterations taken to converge.

    Notes
    -----
    See examples/cluster/plot_affinity_propagation.py for an example.

    The algorithmic complexity of affinity propagation is quadratic
    in the number of points.

    References
    ----------
    """

    def __init__(self, damping=.5, max_iter=200, convergence_iter=15,
                 copy=True, preference=None, affinity=np.corrcoef,
                 verbose=False):

        self.damping = damping
        self.max_iter = max_iter
        self.convergence_iter = convergence_iter
        self.copy = copy
        self.verbose = verbose
        self.preference = preference
        self.affinity = affinity

    def fit(self, X, max_cutoff, min_cutoff):
        """ Create affinity matrix from negative euclidean distances, then
        apply affinity propagation clustering.

        Parameters
        ----------

        X: array-like, shape (n_samples, n_features) or (n_samples, n_samples)
            Data matrix or, if affinity is ``precomputed``, matrix of
            similarities / affinities.
        """



        # TO DO
        X = np.asarray(X)

        self.data = X

        self.affinity_matrix_ = self.affinity(X)

        ####

        self.cluster_centers_indices_, self.labels_, self.n_iter_ = \
            affinity_propagation(X,
                self.affinity_matrix_, max_cutoff, min_cutoff,
                self.preference, max_iter=self.max_iter, convergence_iter=self.convergence_iter,
                damping=self.damping, copy=self.copy, verbose=self.verbose, return_n_iter=True)

        return self

    def predict(self, X):
        """Predict the closest cluster each sample in X belongs to.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            New data to predict.

        Returns
        -------
        labels : array, shape (n_samples,)
            Index of the cluster each sample belongs to.
        """
        ## TO DO
        pass

if __name__ == "__main__":
    cluster = DistinctAffinityPropagation()
    df = pd.read_csv("./csv/chr3_187450000_187470000.csv", sep="\t")
    data = df.as_matrix()
    sample_names = data[:, 0]
    data_values = data[:, 1:].astype(np.float64)

    cluster.fit(data_values, 0.5, 0)

    np.set_printoptions(threshold=np.nan)

    # print cluster.affinity_matrix_

    print cluster.cluster_centers_indices_