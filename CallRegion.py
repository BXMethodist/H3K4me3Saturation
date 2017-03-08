# This is the main module to initiate the region calling

from Wig import Wig
import pickle, numpy as np, pandas as pd
from collections import defaultdict
from multiprocessing import Process, Queue
from predict import optimize_allocs
from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit, Annotation


def CallRegion(wigs, refmap, genome_size_path, process=8):
    """
    :param wig: a dictionary contains wig file paths. key is group name, value is the list of wig file paths in the group
    :param refmap: the path for the reference map
    :param callvariant: whether to call variants expression level
    :return:
    """

    # Load refmap
    with open(refmap, 'rb') as f:
        region_map = pickle.load(f)
    f.close()

    dfs_region = []
    dfs_variant =[]
    for key, value in wigs.items():
        for cur_wig in value:
            region, variant = CallVariants(cur_wig, region_map, process)
            df_region = pd.DataFrame(region)
            df_variant = pd.DataFrame(variant)
            dfs_region.append(df_region)
            dfs_variant.append(df_variant)
            df_region.to_csv(key+'_region.csv', header=False, index=False)
            df_variant.to_csv(key + '_variant.csv', header=False, index=False)

def CallVariants(wig, refmap, process):
    """
    :param wig: the Wig object
    :param refmap: referencemap path, a dictionary group by chromosome
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

    region_results = []
    variant_results = []

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
        cur_region_result, cur_variant_result = queue.get()
        region_results += cur_region_result
        variant_results += cur_variant_result

    for p in processes:
        p.join()

    return region_results, variant_results

def CallVariantsProcess(wigchrome, refmap, queue):
    cur_region_results = []
    cur_variant_results = []
    for key in refmap.keys():
        cur_chrmap = refmap[key]
        cur_wigchrome = wigchrome[key]
        for region in cur_chrmap:
            cur_data = cur_wigchrome.get_signals(region.start, region.end)
            if region.plot:
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

                cur_region_results.append((region.id, cur_total_signals, predict_signals, error))
            else:
                cur_total_signals = np.sum(cur_data)
                predict_signals = cur_total_signals
                error = abs(cur_total_signals - predict_signals) / cur_total_signals
                cur_region_results.append((region.id, cur_total_signals, predict_signals, error))
    queue.put((cur_region_results, cur_variant_results))
    return


Annotation('./75refmap_combined_3kb_regions.pkl','75_combined_3kb')

with open('./superwig.pkl', 'rb') as f:
    superwig = pickle.load(f)
f.close()

wigs = {'super':[superwig]}

path = './75_combined_3kb.pkl'
genomesize = '/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt'

CallRegion(wigs, path, genomesize)