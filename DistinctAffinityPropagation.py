import numpy as np, pandas as pd, scipy
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity
from clusterUtils import *


def affinity_propagation(X, S, max_cutoff, min_cutoff, number_sample_used,
                         preference=None, convergence_iter=15, max_iter=200,
                         damping=0.5, copy=True, verbose=False,
                         return_n_iter=False):
    """Perform Affinity Propagation Clustering of data

    Read more in the :ref:`User Guide <affinity_propagation>`.

    Parameters
    ----------

    X : array-like, shape (n_samples, n_features)
        Original samples matrix

    S : array-like, shape (n_samples, n_samples)
        Matrix of similarities between points

    max_cutoff: the min cut_off to group sample together by affinity metrics

    min_cutoff: the max cut_off to find the most distinct pair by affinity metrics

    number_sample_used: total sample used to build the reference map

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

    seeds : array, shape (n_clusters, n_features)
        the final representation of the clusters

    Notes
    -----
    This function should only be called by .fit function

    References
    ----------
    """

    n_samples = S.shape[0]
    cluster_centers_indices = []
    labels = []
    seeds = []
    frontiers = np.zeros(n_samples)
    n_iter = 0

    affinity_matrix = S.copy()

    if S.shape[0] != S.shape[1]:
        raise ValueError("S must be a square array (shape=%s)" % repr(S.shape))

    while len(cluster_centers_indices) == 0:
        distinct_index = np.argmin(affinity_matrix)

        min_distance = affinity_matrix.flat[distinct_index]

        if min_distance > min_cutoff or np.isnan(min_distance) or S.shape[0] < number_sample_used/5:
            # TO DO: return the entire samples as one cluster
            labels.append(np.arange(n_samples))
            representation = [np.mean(X, axis=0)]
            seeds = [np.mean(X, axis=0)]
            cluster_centers_indices = [0]
            return cluster_centers_indices, labels, seeds, representation

        distinct_pairs = np.asarray(np.where(affinity_matrix == min_distance)).T

        distinct_set = remove_duplicate(distinct_pairs)

        if len(distinct_set) == 1:
            first_distinct = distinct_set[0]
        else:
            first_distinct = most_different_pair(X, distinct_set)

        # print first_distinct
        for i in range(len(first_distinct)):
            affinity_x = affinity_matrix[first_distinct[i], :]
            cluster_x = np.intersect1d(np.where(affinity_x > max_cutoff), np.where(affinity_x <= 1))
            if len(cluster_x) >= 5 and len(cluster_x) >= number_sample_used/50:
                labels.append(cluster_x)
                cluster_centers_indices.append(first_distinct[i])
                seeds.append(X[first_distinct[i]])

                frontiers[cluster_x] = 1
                affinity_matrix[cluster_x, :] = float("inf")
                affinity_matrix[:, cluster_x] = float("inf")
            else:
                frontiers[first_distinct[i]] = 1
                affinity_matrix[first_distinct[i], :] = float("inf")
                affinity_matrix[:, first_distinct[i]] = float("inf")

    while True and n_iter < max_iter:
        samples_left = np.where(frontiers==0)[0]

        new_centers = []

        for sample in samples_left:
            distinct_qualifier = [S[sample, x] <= min_cutoff for x in cluster_centers_indices]

            if all(distinct_qualifier):
                new_centers.append(sample)

        if len(new_centers) == 0:
            break
        elif len(new_centers) > 1:
            # TO DO: create new function to get better performance
            best_new_center = accumulate_selecter(X, cluster_centers_indices, new_centers, S)
        else:
            best_new_center = new_centers[0]

        affinity_new_center = affinity_matrix[best_new_center, :]
        cluster_new_center = np.intersect1d(np.where(affinity_new_center > max_cutoff),
                                            np.where(affinity_new_center <= 1))

        if len(cluster_new_center) >= 5 and len(cluster_new_center) >= S.shape[0]/50:
            labels.append(cluster_new_center)
            cluster_centers_indices.append(best_new_center)
            seeds.append(X[best_new_center])

            frontiers[cluster_new_center] = 1

            affinity_matrix[cluster_new_center, :] = float("inf")
            affinity_matrix[:, cluster_new_center] = float("inf")

        else:
            frontiers[best_new_center] = 1
            affinity_matrix[best_new_center, :] = float("inf")
            affinity_matrix[:, best_new_center] = float("inf")

        if np.sum(frontiers) == n_samples:
            break
        n_iter += 1

    if len(labels) == 0:
        labels.append(np.arange(n_samples))
        representation = [np.mean(X, axis=0)]
        seeds = [np.mean(X, axis=0)]
        cluster_centers_indices = [0]
        return cluster_centers_indices, labels, seeds, representation

    representation = [np.mean(X[label], axis=0) for label in labels]
    representation = np.asarray(representation)

    # enriched_representation = []
    # for i in range(representation.shape[0]):
    #     max_value = np.amax(representation[i, :])
    #     enriched_rep = representation[i, :] * representation[i, :] / max_value
    #     enriched_representation.append(enriched_rep)
    #
    # feature_mark = np.argmax(enriched_representation, axis=0)
    # for i in range(len(representation)):
    #     representation[i][feature_mark!=i] = 0


    ### use representation to recluster
    reclusters = labels

    # for i in range(representation.shape[0]):
    #     new_seed = representation[i, :]
    #     similarity = cosine_similarity(X, new_seed).T[0]
    #     reclusters[i] = np.where(similarity > max_cutoff)[0]
    representation = mean_normalization(representation)

    return cluster_centers_indices, reclusters, seeds, representation

    # return cluster_centers_indices, labels, n_iter

def signal_reduction(sample, array1, array2):
    signal_enrich = np.zeros(sample.shape[0])
    result1 = sample - array1
    result1[result1<0] = 0
    signal_enrich += np.sum(result1, axis=0)
    result2 = sample - array2
    result2[result2 < 0] = 0
    signal_enrich += np.sum(result2, axis=0)

    return signal_enrich

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

def accumulate_selecter(samples, cluster_centers_indices, candidates, affinity_matrix):
    """Make sure the new cluster are also very different with the old ones.

                Parameters
                ----------
                samples: ndarray, (nsample, nfeatures)
                cluster_centers_indices: existing clusters index
                candidates: potential new cluster index
                affinity_matrix: distance matrix for all the samples, original one, not modified.
                Return
                ----------
                the most distinct candidates

                Notes
                -----

                References
                ----------
                """
    best_candidates = []
    min_dis = None
    for candidate in candidates:
        cur_best = np.max(affinity_matrix[candidate, cluster_centers_indices])
        if min_dis is None:
            min_dis = cur_best
            best_candidates.append(candidate)
        elif min_dis == cur_best:
            best_candidates.append(candidate)
        elif min_dis > cur_best:
            min_dis = cur_best
            best_candidates = [candidate]
    if len(best_candidates) == 1:
        return best_candidates[0]

    candidates = best_candidates

    best_candidates_total_dis = []
    total_dis = None
    for candidate in candidates:
        cur_best = np.sum(affinity_matrix[candidate, cluster_centers_indices])
        if total_dis is None:
            total_dis = cur_best
            best_candidates_total_dis.append(candidate)
        elif total_dis == cur_best:
            best_candidates_total_dis.append(candidate)
        elif total_dis > cur_best:
            total_dis = cur_best
            best_candidates_total_dis = [candidate]

    if len(best_candidates_total_dis) == 1:
        return best_candidates_total_dis[0]

    candidates = best_candidates_total_dis

    max_distance = float("-inf")
    best_candidate = None
    for candidate in candidates:
        cur_distance = 0
        for center in cluster_centers_indices:
            cur_distance += euclidean(samples[candidate], samples[center])
        if cur_distance > max_distance:
            best_candidate = candidate
            max_distance = cur_distance
    return best_candidate


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

    def __init__(self, number_sample_used, damping=.5, max_iter=200, convergence_iter=15,
                 copy=True, preference=None, affinity=np.corrcoef,
                 verbose=False):

        self.damping = damping
        self.max_iter = max_iter
        self.convergence_iter = convergence_iter
        self.copy = copy
        self.verbose = verbose
        self.preference = preference
        self.affinity = affinity
        self.number_sample_used = number_sample_used

    def fit(self, X, max_cutoff, min_cutoff, mean_target=200, smooth_period=20, cutoff=None):
        """ Create affinity matrix from negative euclidean distances, then
        apply affinity propagation clustering.

        Parameters
        ----------

        X: array-like, shape (n_samples, n_features) or (n_samples, n_samples)
            Data matrix or, if affinity is ``precomputed``, matrix of
            similarities / affinities.
        """
        X = np.asarray(X)
        # print X
        # if mean_target:
        #     X = mean_normalization(X, target_signal=mean_target)
        if smooth_period:
            X = smooth_normalization(X, smooth=smooth_period)

        self.data = X

        self.affinity_matrix_ = self.affinity(X)

        ####

        self.cluster_centers_indices_, self.labels_, self.seeds_, self.representation_ = \
            affinity_propagation(X,
                self.affinity_matrix_, max_cutoff, min_cutoff, self.number_sample_used,
                self.preference, max_iter=self.max_iter, convergence_iter=self.convergence_iter,
                damping=self.damping, copy=self.copy, verbose=self.verbose, return_n_iter=True)

        return self

