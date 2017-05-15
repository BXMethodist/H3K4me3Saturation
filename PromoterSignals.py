"""
This module is to check the signals around promoter, to see why they are not covered by the reference map
"""

from Wig import Wig
import pickle, os, pandas as pd, numpy as np
from refMapUtils import load_obj
from wigChrom import WigChrom
from collections import defaultdict
from multiprocessing import Process, Queue

def gene_tss_signals(distance, gene_distance_table, wigs_obj_path, gene_gtf, process=40):
    wigs_obj_path += '' if wigs_obj_path.endswith('/') else '/'
    wigs = [wigs_obj_path+x for x in os.listdir(wigs_obj_path)]

    chunks = []
    cur_index = 0
    reminder = len(wigs) % process
    chunk_size = len(wigs) / process
    for i in range(process):
        if reminder > 0:
            chunks.append(wigs[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
            cur_index += 1
            reminder -= 1
        else:
            chunks.append(wigs[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

    total_chunk_size = 0
    for chunk in chunks:
        total_chunk_size += len(chunk)
    if total_chunk_size != len(wigs):
        print 'multiple processes chunk size is not correct'
        return None

    queue = Queue()
    processes = []

    for i in range(process):
        cur_chunk = chunks[i]
        p = Process(target=gene_tss_signal,
                    args=(queue, distance, gene_distance_table, cur_chunk, gene_gtf))
        processes.append(p)
        p.start()

    gene_promoter_signals = {}

    for i in range(process):
        cur_signals = queue.get()
        for key, value in cur_signals.items():
            if key not in gene_promoter_signals:
                gene_promoter_signals[key] = value
            else:
                prev = gene_promoter_signals[key]
                new_sig = (prev[1]+value[1])/2.0
                new_height = (prev[2] + value[2])/2.0
                gene_promoter_signals[key]=[prev[0], new_sig, new_height, prev[-1]]
    for p in processes:
        p.join()

    results = []
    for key, value in gene_promoter_signals.items():
        results.append(value)

    result_df = pd.DataFrame(results)
    result_df.columns=['gene_name2', 'signals', 'heights', 'group']
    result_df.set_index(['gene_name2'])

    result_df.to_csv('gene_promoter_signals_100_10.csv')

    return result_df

def gene_tss_signal(queue, distance, gene_distance_table, wigs, gene_gtf):
    """

    :param distance: distance cutoff to check the signals
    :param gene_distance_table: gene distance table towards reference map, 3rd column is the total distance
    :param wigs_obj_path: wigs object path
    :return:
    """

    df = pd.read_csv(gene_distance_table, index_col=0, header=None)
    df.columns = ['up','down','total']

    gtf_df = pd.read_csv(gene_gtf, sep='\t')

    gene_signals = {}

    for wig in wigs:
        cur_wig = load_obj(wig[:-4])

        for gene_name2 in df.index:
            cur_df = gtf_df[gtf_df['name2'] == gene_name2]
            chromosome = cur_df.iloc[0]['chrom']
            best_signal, best_height = 0, 0

            if chromosome not in cur_wig.genome:
                print gene_name2
                continue

            for i in range(cur_df.shape[0]):
                if cur_df.iloc[i]['strand'] == '+':
                    cur_tss = int(cur_df.iloc[i]['txStart'])
                elif cur_df.iloc[i]['strand'] == '-':
                    cur_tss = int(cur_df.iloc[i]['txEnd'])
                start, end = cur_tss-distance, cur_tss+distance
                cur_signals_vecter = cur_wig.genome[chromosome].get_signals(start, end)

                best_signal = max(np.sum(cur_signals_vecter), best_signal)
                best_height = max(np.max(cur_signals_vecter), best_height)

            if df.ix[gene_name2, 'total'] > distance:
                cur_group ='uncovered'
            else:
                cur_group = 'covered'
            if gene_name2 not in gene_signals:
                gene_signals[gene_name2] = [gene_name2, best_signal, best_height, cur_group]
            else:
                prev = gene_signals[gene_name2]
                new_signals = prev[1] + best_signal
                new_heights = prev[2] + best_height
                gene_signals[gene_name2] = [gene_name2, new_signals, new_heights, cur_group]

    result_gene_signals = {}
    sample_numbers = len(wigs)

    for key, value in gene_signals.items():
        cur_signals, cur_heights = value[1], value[2]
        cur_signals = cur_signals*1.0/sample_numbers
        cur_heights = cur_heights*1.0/sample_numbers
        result_gene_signals[key] = ([value[0], cur_signals, cur_heights, value[-1]])

    queue.put(result_gene_signals)
    return

# gene_tss_signals(10000, './100_extend_10_map/100_10_gene_distance.csv', './wigs', './pkl/hg19_RefSeq_refGene.txt')

