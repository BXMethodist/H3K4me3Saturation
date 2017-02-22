import scipy.optimize as spo
from scipy.spatial.distance import euclidean
import numpy as np


def optimize_allocs(data, clusters):
    """
    :param data: ndarray for signal need to be predited
    :param clusters: clusters from refmap
    :return: allocs, a list of quotes of clusters
    """

    allocs = np.asarray([1.0/clusters.shape[0]]*clusters.shape[0])

    allocs = spo.minimize(distance, allocs, args=(data, clusters), method='SLSQP').x

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
    return dist

