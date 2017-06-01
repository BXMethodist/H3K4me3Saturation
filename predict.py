import scipy.optimize as spo
from scipy.spatial.distance import euclidean
import numpy as np


def optimize_allocs(data, clusters):
    """
    :param data: ndarray for signal need to be predited
    :param clusters: clusters from refmap
    :return: allocs, a list of quotes of clusters
    """
    clusters = np.asarray(clusters)

    norm_factor = np.sum(data)

    new_clusters = []

    for i in range(clusters.shape[0]):
        c = clusters[i, :]
        c = c * norm_factor/np.sum(c)
        new_clusters.append(c)
    clusters = np.asarray(new_clusters)

    constraints = {'type': 'ineq', 'fun': lambda x: x - 0.0}

    total_signals = np.sum(data)

    variants_signals = np.sum(clusters, axis=1)

    borders = total_signals/variants_signals

    bounds = tuple([(0.0, 1.0) for i in range(len(borders))])

    # options = {'disp': False, 'maxiter': 10000}
    options = {'disp': False}
    allocs = np.asarray([1.0/clusters.shape[0]]*clusters.shape[0])

    allocs = spo.minimize(distance, allocs, args=(data, clusters),
                          method='SLSQP', constraints=constraints, bounds=bounds, options=options).x

    # print allocs
    # allocs = allocs/np.sum(allocs)
    return allocs


def distance(allocs, data, clusters):
    """
    :param data: ndarray contains data to be analyzed.
    :param clusters: clusters in refmap, np array
    :param allocs: multiplication factors for each clusters
    :return: distance towards target
    """
    predict_signal = np.zeros(clusters.shape[1])
    for i in range(len(allocs)):
        quote = allocs[i]
        cluster = clusters[i]
        predict_signal += quote*cluster
    dist = euclidean(data, predict_signal)
    # dist = abs(np.corrcoef(data, predict_signal)[0, 1] - 1.0)
    return dist

