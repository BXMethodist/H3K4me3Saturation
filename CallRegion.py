# This is the main module to initiate the region calling

from predict import optimize_allocs
from Wig import Wig
import pickle, numpy as np, pandas as pd


def CallRegion(wigs, refmap, genome_size_path, callvariant=True, process=8):
    """
    :param wig: a dictionary contains wig file paths. key is group name, value is the list of wig file paths in the group
    :param refmap: the path for the reference map
    :param callvariant: whether to call variants expression level
    :return:
    """

    dfs = []
    for key, value in wigs.items():
        for path in value:
            cur_wig = Wig(path, genome_size_path)
            df = CallVariants(cur_wig, refmap, key, callvariant, process)
            dfs.append(df)


def CallVariants(wig, refmap, group, callvariant, process):
    """
    :param wig: the Wig object
    :param refmap: referencemap path, a dictionary group by chromosome
    :param group: the group of sample belongs to
    :param callvariant: whether to call variant
    :param process: number of process
    :return: dataframe
    """
    with open(refmap, 'rb') as f:
        reference_map = pickle.load(f)

    results = []
    genome_signals = wig.genome
    if callvariant:
        for region in reference_map:
            chromosome, start, end = region.chromosome, region.start, region.end
            target = genome_signals[chromosome].get_signals(start, end)
            allocs = optimize_allocs(target, region.representatives)
            total_target_signals = np.sum(target)
            allocs = np.asarray(allocs)/np.sum(allocs)
            cur_result = []
            for i in range(len(region.representatives)):
                name = "_".join([chromosome, str(start), str(end), 'variant '+str(i)])
                cur_result.append([name, total_target_signals*allocs[i]])
            results += cur_result
    else:
        for region in reference_map:
            chromosome, start, end = region.chromosome, region.start, region.end
            target = genome_signals[chromosome].get_signals(start, end)
            total_target_signals = np.sum(target)
            cur_result = []
            for i in range(len(region.representatives)):
                name = "_".join([chromosome, str(start), str(end), 'variant '+str(i)])
                cur_result.append([name, total_target_signals])
            results += cur_result

    df = pd.DataFrame(results, columns=['variant name', 'value'], index='variant')
    return df

def CallUnits():
    pass
