# This is the main module to initiate the region calling

from Wig import Wig
import pickle, numpy as np, pandas as pd, os
from collections import defaultdict
from multiprocessing import Process, Queue
from predict import optimize_allocs
from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit, Annotation
import scipy.stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

def CallRegion(wigs, refmap, genome_size_path, output, alias=None, process=8):
    """
    :param wig: a dictionary contains wig file paths. key is group name, value is the list of wig file paths in the group
    :param refmap: the path for the reference map
    :param callvariant: whether to call variants expression level
    :param alias: a map contains group - samples alias for output table
    :return: variant df and region df
    """

    # Load refmap
    with open(refmap, 'rb') as f:
        region_map = pickle.load(f)
    f.close()

    # for key, value in wigs.items():
    #     new_wig_objs = []
    #     for wig in value:
    #         cur_wig = Wig(wig, genome_size_path)
    #         new_wig_objs.append(cur_wig)
    #     wigs[key] = new_wig_objs

    rownames_region = [region.id for region in region_map]
    rownames_variant = [variant.id for region in region_map for variant in region.variants]

    dfs_region_error = pd.DataFrame(index=rownames_region)
    dfs_variant =pd.DataFrame(index=rownames_variant)
    dfs_region = pd.DataFrame(index=rownames_region)

    groupnames = defaultdict(list)

    for key, value in wigs.items():
        for i in range(len(value)):
            cur_wig = value[i]
            colname = key+'_'+cur_wig.file_name if alias is None else alias[key][i]
            groupnames[key].append(colname)
            region, region_error, variant = CallVariants(cur_wig, region_map, process)

            df_region_error = pd.DataFrame(region_error,
                                     columns=['region_id', colname+'_target', colname+"_predict", colname+"_error"])
            df_region_error = df_region_error.set_index(['region_id'])

            df_region = pd.DataFrame(region, columns=['region_id', colname])
            df_region = df_region.set_index(['region_id'])

            df_variant = pd.DataFrame(variant,
                                      columns=['variant_id', colname])
            df_variant = df_variant.set_index(['variant_id'])

            dfs_region_error = dfs_region_error.join(df_region_error)
            dfs_variant = dfs_variant.join(df_variant)
            dfs_region = dfs_region.join(df_region)

    for key in groupnames.keys():
        dfs_variant[key] = dfs_variant[groupnames[key]].mean(axis=1)
        dfs_region[key] = dfs_region[groupnames[key]].mean(axis=1)

    min_variant = min(dfs_variant[[key for key in groupnames.keys()]].min())
    min_region = min(dfs_region[[key for key in groupnames.keys()]].min())

    for key in groupnames.keys():
        dfs_variant[key] = dfs_variant[key] + min_variant
        dfs_region[key] = dfs_region[key] + min_region

    stats = importr('stats')
    for i in range(len(groupnames.keys())):
        key1 = groupnames.keys()[i]
        for j in range(i+1, len(groupnames.keys())):
            key2 = groupnames.keys()[j]
            dfs_variant[key1+"_vs_"+key2+"_log2FC"] = np.log2(dfs_variant[key1]/dfs_variant[key2])
            dfs_variant[key1 + '_vs_' + key2 + "_P"] = 1 - scipy.stats.poisson.cdf(
                dfs_variant[[key1, key2]].max(axis=1),
                dfs_variant[[key1, key2]].min(axis=1))
            dfs_variant[key1 + '_vs_' + key2 + "_log10P"] = np.log10(dfs_variant[key1 + '_vs_' + key2 + "_P"])
            dfs_variant[key1 + '_vs_' + key2 + "_FDR"] = stats.p_adjust(
                FloatVector(dfs_variant[key1 + '_vs_' + key2 + "_P"].tolist()),
                method='BH')


            dfs_region[key1 + "_vs_" + key2 + "_log2FC"] = np.log2(dfs_region[key1] / dfs_region[key2])
            dfs_region[key1 + '_vs_' + key2 + "_P"] = 1-scipy.stats.poisson.cdf(
                                                                dfs_region[[key1, key2]].max(axis=1),
                                                                dfs_region[[key1, key2]].min(axis=1))
            dfs_region[key1 + '_vs_' + key2 + "_log10P"] = np.log10(dfs_region[key1 + '_vs_' + key2 + "_P"])
            dfs_region[key1 + '_vs_' + key2 + "_FDR"] = stats.p_adjust(
                FloatVector(dfs_region[key1 + '_vs_' + key2 + "_P"].tolist()),
                method='BH')

    dfs_region_error.to_csv(output+'_region_error.csv')
    dfs_variant.to_csv(output + '_variant.csv')
    dfs_region.to_csv(output + '_region.csv')
    return dfs_variant, dfs_region


def CallVariants(wig, refmap, process):
    """
    :param wig: the Wig object
    :param refmap: referencemap, a dictionary group by chromosome
    :param group: the group of sample belongs to
    :param callvariant: whether to call variant
    :param process: number of process
    :return: dataframe
    """
    chromosomes = wig.genome.keys()
    chunk_size = len(chromosomes)/process
    reminder = len(chromosomes)%process

    regionmap = defaultdict(list)
    for region in refmap:
        regionmap[region.chromosome].append(region)

    chunks = []
    cur_index = 0
    for i in range(process):
        if reminder > 0:
            chunks.append(chromosomes[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
            cur_index += 1
            reminder -= 1
        else:
            chunks.append(chromosomes[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

    queue = Queue()
    processes = []

    region_error_results = []
    variant_results = []
    region_results = []

    # print chunks

    for i in range(process):
        cur_chrs = chunks[i]
        cur_refmap = {}
        cur_wig = {}
        for cur_chr in cur_chrs:
            cur_refmap[cur_chr] = regionmap[cur_chr]
            cur_wig[cur_chr] = wig.genome[cur_chr]
        p = Process(target=CallVariantsProcess, args=(cur_wig, cur_refmap, queue))
        processes.append(p)
        p.start()

    for i in range(process):
        cur_region, cur_region_result_error, cur_variant_result = queue.get()
        region_results += cur_region
        region_error_results += cur_region_result_error
        variant_results += cur_variant_result

    for p in processes:
        p.join()

    return region_results, region_error_results, variant_results

def CallVariantsProcess(wigchrome, refmap, queue):
    # print os.getpid()
    cur_region_results = []
    cur_region_results_error = []
    cur_variant_results = []
    for key in refmap.keys():
        cur_chrmap = refmap[key]
        cur_wigchrome = wigchrome[key]
        # n = 0
        for region in cur_chrmap:
            cur_data = cur_wigchrome.get_signals(region.start, region.end)
            # print "fetch data complete for ", region.id
            # print n
            # n +=1

            cur_variant_representatives = []
            cur_ids = []
            for variant in region.variants:
                cur_ids.append(variant.id)
                cur_variant_representatives.append(variant.representative)
            cur_allocs = optimize_allocs(cur_data, cur_variant_representatives)

            cur_total_signals = np.sum(cur_data)
            predict_signals = 0

            for i in range(len(region.variants)):
                cur_var_signal = cur_total_signals*cur_allocs[i]
                cur_variant_results.append((cur_ids[i], cur_var_signal))
                predict_signals += cur_var_signal
            error = abs(cur_total_signals-predict_signals)/cur_total_signals
            # print (region.id, cur_total_signals, predict_signals, error)
            cur_region_results_error.append((region.id, cur_total_signals, predict_signals, error))
            cur_region_results.append((region.id, cur_total_signals))

    queue.put((cur_region_results, cur_region_results_error, cur_variant_results))
    return

def DiffVariant(refmap, dfs_variant, dfs_region):
    """
    this function is used to call the pattern alterations
    :param refmap: reference map, a dictionary group by chromosome
    :param dfs_variant: data frame from callregion for variant
    :param dfs_region: data frame from callregion for region
    :return: a df containing all the variant and regions that involved in pattern alteration.
    """
    # to do:
    pass

# Annotation('./75refmap_combined_3kb_regions.pkl','75_combined_3kb')

with open('./wig/superwig.pkl', 'rb') as f:
    superwig = pickle.load(f)
f.close()
print "loading complete"
wigs = {'super1':[superwig], 'super2':[superwig]}
# print superwig.genome['chr4'].get_signals(9980, 10280)

path = './pkl/75_combined_3kb.pkl'
genomesize = '/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt'
#
CallRegion(wigs, path, genomesize, 'super', process=8)