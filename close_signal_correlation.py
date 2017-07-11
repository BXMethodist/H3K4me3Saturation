import os, pandas as pd, numpy as np, pickle, random, bisect, random
from refMapUtils import load_obj
from collections import defaultdict


def closest_peak(seed, peaks):
    indexes = [p[1] for p in peaks]
    index = bisect.bisect(indexes, seed[1])

    if peaks[index-1][1] != seed[1]:
        print peaks[index-1],peaks[index], seed
        print "something wrong?!"
    if index-2 >= 0:
        left = peaks[index-2]
        # print left, seed
        left_distance = seed[1] - left[2]
    else:
        left_distance = None
        left = None

    if index < len(peaks):
        right = peaks[index]
        right_distance = right[1] - seed[2]
    else:
        right_distance = None
        right = None

    if left_distance is None:
        return right
    elif right_distance is None:
        return left
    elif left_distance < right_distance:
        return left
    elif right_distance <= left_distance:
        return right

def get_wig_signals(wig, regions):
    cur_wig = load_obj(wig[:-4])
    results = []
    for region in regions:
        if region[0] in cur_wig.genome:
            cur_signals = cur_wig.genome[region[0]].get_signals(region[1], region[2])
            # print region
            # print cur_signals
            # print type(region[0]), type(region[1])
            results.append(float(format(np.sum(cur_signals), '.2f')))
        else:
            results.append(0.0)
    return results

def close_peak_signals(peak_tables):
    sample_peak_numbers = {}
    total_number = 0
    sample_peaks = {}
    sample_peaks_list = {}

    all_peaks = defaultdict(set)

    chr_df = pd.read_csv('./pkl/hg19_chr_sizes.txt', sep='\t', index_col=0, header=None)
    chromosomes = [chr.strip() for chr in chr_df.index]

    for i in range(len(peak_tables)):
        peak_table = open(peak_tables[i], 'r')
        info1 = peak_table.readlines()
        peak_table.close()

        peaks = {}
        peaks_list = set()

        for line1 in info1:
            line1 = line1.strip()
            line1 = line1.split('\t')
            if line1[1] == "start":
                continue

            chromosome = line1[0]
            if chromosome not in chromosomes:
                continue

            if chromosome not in peaks:
                peaks[chromosome] = set()
            peaks[chromosome].add((chromosome, int(line1[1]), int(line1[2])))
            peaks_list.add((chromosome, int(line1[1]), int(line1[2])))

        for chromosome in peaks.keys():
            all_peaks[chromosome] = all_peaks[chromosome].union(peaks[chromosome])

        sample_peaks[i] = peaks
        total_number += len(peaks_list)
        sample_peak_numbers[i] = len(peaks_list)
        sample_peaks_list[i] = list(peaks_list)

    random_size = 100000

    correlations_regions = set()
    correlations_pairs = set()

    for chromosome in all_peaks.keys():
        all_peaks[chromosome] = sorted(list(all_peaks[chromosome]), key=lambda x:x[1])

    for i in range(len(peak_tables)):
        cur_random_size = int(sample_peak_numbers[i]*1.0/total_number * random_size)
        cur_random_peaks_indexes = random.sample(range(0, len(sample_peaks_list[i])), cur_random_size)

        cur_random_seeds = [sample_peaks_list[i][j] for j in cur_random_peaks_indexes]

        for cur_seed in cur_random_seeds:
            cur_closet = closest_peak(cur_seed, all_peaks[cur_seed[0]])
            if cur_closet is not None:
                correlations_regions.add(cur_seed)
                correlations_regions.add(cur_closet)
                if ('_'.join(list([str(x) for x in cur_seed])), '_'.join(list([str(y) for y in cur_closet]))) not in correlations_pairs \
                    and ('_'.join(list([str(x) for x in cur_closet])), '_'.join(list([str(y) for y in cur_seed]))) not in correlations_pairs:
                    correlations_pairs.add(('_'.join(list([str(x) for x in cur_seed])), '_'.join(list([str(y) for y in cur_closet]))))
            else:
                continue

    wigs_signals = []
    correlations_regions = list(correlations_regions)
    correlations_pairs = list(correlations_pairs)

    wigs = ['./wigs/'+ x for x in os.listdir('./wigs/')]

    # correlations_regions = [('chr4', 44728720, 44728960)]

    for wig in wigs:
        wigs_signals.append(get_wig_signals(wig, correlations_regions))


    df = pd.DataFrame(wigs_signals)
    df.columns= ['_'.join(list([str(z) for z in region])) for region in correlations_regions]

    df.to_csv('peaks_correlation.csv', index=None)

    df_cols = pd.DataFrame(correlations_pairs)

    df_cols.to_csv('peaks_correlation_columns.csv', index=None, header=None)

    return df, df_cols

def read_refmap(refmap, step=10):
    """
    read refmap csv file, change it to a hashmap, with real genome index
    :param refmap:
    :return:
    """
    f = open(refmap, 'r')
    regions = []
    regions_dict = defaultdict(list)
    cur_chr = None
    for line in f:
        if line.startswith(">"):
            cur_chr = line.strip()[1:]
        elif len(line.strip()) > 1:
            regions.append([cur_chr] + [int(x) * step for x in line.strip().split(',')])
            regions_dict[cur_chr].append([cur_chr] + [int(x) * step for x in line.strip().split(',')])
    return regions, regions_dict

def region_signals(refmap_list, wigs):
    wigs_signals = []

    for wig in wigs:
        wigs_signals.append(get_wig_signals(wig, refmap_list))

    df = pd.DataFrame(wigs_signals)
    df.columns = ['_'.join(list([str(z) for z in region])) for region in refmap_list]

    df.to_csv('region_correlation.csv', index=None)
    return df

def region_random_pair(refmap_dict, random_number=1000000):
    total_peaks = 0
    region_count = defaultdict(int)
    for key, value in refmap_dict.items():
        total_peaks += len(value)
        region_count[key] = len(value)
    for key, value in region_count.items():
        qouta = value * 1.0 / total_peaks * random_number
        if qouta <= value:
            region_count[key] = qouta
        else:
            pass
    region_pair = []
    for key, value in refmap_dict.items():
        for i in range(region_count[key]):
            indices = random.sample(range(len(value)), 2)

            p1 = '_'.join([str(x) for x in value[indices[0]]])
            p2 = '_'.join([str(x) for x in value[indices[1]]])

            region_pair.append((key, p1, p2))
    df = pd.DataFrame(region_pair)
    df.to_csv('region_correlation_pair.csv', index=None, header=None)
    return df

def peaks_correlation(signals_table, cols_table, random_pairs):
    cols_df = pd.read_csv(cols_table, header=None)
    print cols_df.shape
    if random_pairs is not None:
        pairs_df = cols_df.sample(random_pairs)
    else:
        pairs_df = cols_df.copy()
    results = defaultdict(set)

    chr_df = pd.read_csv('./pkl/hg19_chr_sizes.txt', sep='\t', index_col=0, header=None)
    chromosomes = [chr.strip() for chr in chr_df.index]
    # print chromosomes

    for i in range(pairs_df.shape[0]):
        p1 = pairs_df.iloc[i, 1]
        p2 = pairs_df.iloc[i, 2]
        p1_str = p1.split('_')
        p2_str = p2.split('_')

        if p1_str[0] not in chr_df.index:
            print p1_str[0]
            continue
        # if p1_str[0].find('chrM') != -1:
        #     print p1_str[0]
        # print p1_str
        if int(p1_str[1]) < int(p2_str[1]):
            distance = int(p2_str[1]) - int(p1_str[2])
        else:
            distance = int(p1_str[1]) - int(p2_str[2])
        # if distance < 0:
        #     print p1, p2
        #     continue
        # else:
        results[distance/100].add((p1, p2, distance))

    signals_df = read_table(signals_table)

    distance_corr_results = defaultdict(list)
    # print results
    print "table load complete!"
    for key in results.keys():
        for p in results[key]:
            r1 = signals_df[p[0]]
            r2 = signals_df[p[1]]
            cur_cor = np.corrcoef(r1, r2)[0, 1]
            # print r1[r1==0].shape[0]
            # print r1[pd.isnull(r1)].shape[0]
            # print r1
            # print cur_cor, p[2]

            distance_corr_results[p[2]/1000*1000].append(cur_cor)

    final_results = []
    for key in distance_corr_results.keys():
        final_results.append((key, np.mean(distance_corr_results[key]), np.std(distance_corr_results[key]), len(distance_corr_results[key])))

    result_df = pd.DataFrame(final_results)
    #
    result_df.to_csv('distance_signal_correlation.csv', index=None, header=None)
    return

def read_table(signal_table):
    f = open(signal_table, 'r')
    columns = f.readline().strip().split(',')
    results = []
    for line in f:
        info = [float(x) if x !='' else 0 for x in line.strip().split(',')]
        results.append(info)

    df = pd.DataFrame(results)
    df.columns = columns
    return df

# from time import time
# start = time()
# read_table('peaks_correlation.csv')
#
# end = time()
# print end - start
# start = time()
# pd.read_csv('peaks_correlation.csv')
#
# end = time()
# print end - start

# close_peak_signals(['/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/' + x for x in
#                     os.listdir('/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/')])

# close_peak_signals(['/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/' + x for x in os.listdir('/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/')])

peaks_correlation('region_correlation.csv', 'region_correlation_pair.csv', None)
# print float(0)
# print float(format(0.00000, '.2f'))
# print float('')

# refmap_list, refmap_dict = read_refmap('100_refmap.csv')
#
# wigs = ['./wigs/'+ x for x in os.listdir('./wigs/')]
#
# region_random_pair(refmap_dict)
# region_signals(refmap_list, wigs)

