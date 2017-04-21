"""
This module is to check the statistics of the reference map, for example, the distribution of width
"""

from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit
import pandas as pd, pickle, numpy as np
from simple_region import SimpleRegion
def distribution_width(referencelistmap, outputname):
    """
    get the distribution of the width of region
    :param referencelistmap:
    :return:
    """
    with open(referencelistmap, 'rb') as f:
        referencelistmap = pickle.load(f)
    f.close()
    results = []
    for region in referencelistmap:
        results.append(region.end - region.start)
    results = sorted(results)

    results = np.asarray(results)
    length = len(results)
    widths = []
    for l in range(100, 100100, 100):
        widths.append(len(results[results<l])*1.0/length)

    df = pd.DataFrame(widths)

    df.to_csv(outputname, index=False, header=None)
    return df

def rank_width(referencelistmap, outputname):
    """
    get the distribution of the width of region
    :param referencelistmap:
    :return:
    """
    with open(referencelistmap, 'rb') as f:
        referencelistmap = pickle.load(f)
    f.close()
    results = []
    for region in referencelistmap:
        results.append([region.end - region.start])
    results = sorted(results, reverse=True)

    results = np.asarray(results)

    df = pd.DataFrame(results, index=np.arange(1, len(results)+1))

    df.to_csv(outputname, index=False, header=None)
    return df


# rank_width("./ref100/100_refregion.pkl", "100_refmap_width_rank.csv")
# rank_width("./extend_100_10/extend_100_10_simple_region.pkl", "100_10_refmap_width_rank.csv")
