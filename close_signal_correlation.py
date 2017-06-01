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

    chr_df = pd.read_csv('./pkl/hg19_chr_sizes.txt', sep='\t', index_col=0, header=None)
    chromosomes = [chr.strip() for chr in chr_df.index]

    for i in range(len(peak_tables)):
        peak_table = open(peak_tables[i], 'r')
        info1 = peak_table.readlines()
        peak_table.close()

        peaks = {}
        peaks_list = []

        for line1 in info1:
            line1 = line1.strip()
            line1 = line1.split('\t')
            if line1[1] == "start":
                continue

            chromosome = line1[0]
            if chromosome not in chromosomes:
                continue

            if chromosome not in peaks:
                peaks[chromosome] = []
            peaks[chromosome].append((chromosome, int(line1[1]), int(line1[2])))
            peaks_list.append((chromosome, int(line1[1]), int(line1[2])))

        number_peaks = 0

        for chromosome in peaks.keys():
            peaks[chromosome] = sorted(peaks[chromosome], key=lambda x: x[1])
            number_peaks += len(peaks[chromosome])

        sample_peaks_list[i] = peaks_list
        sample_peaks[i] = peaks
        sample_peak_numbers[i] = number_peaks
        total_number += number_peaks

    random_size = 100000

    correlations_regions = set()
    correlations_pairs = set()

    for i in range(len(peak_tables)):
        cur_random_size = int(sample_peak_numbers[i]*1.0/total_number * random_size)
        cur_random_peaks_indexes = random.sample(range(0, len(sample_peaks_list[i])), cur_random_size)

        cur_random_seeds = [sample_peaks_list[i][j] for j in cur_random_peaks_indexes]

        for cur_seed in cur_random_seeds:
            cur_closet = closest_peak(cur_seed, sample_peaks[i][cur_seed[0]])
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

    for wig in wigs[10:20]:
        wigs_signals.append(get_wig_signals(wig, correlations_regions))


    df = pd.DataFrame(wigs_signals)
    df.columns= ['_'.join(list([str(z) for z in region])) for region in correlations_regions]

    df.to_csv('peaks_correlation.csv')

    df_cols = pd.DataFrame(correlations_pairs)

    df_cols.to_csv('peaks_correlation_columns.csv')

    return df, df_cols

def peaks_correlation(signals_table, cols_table, random_pairs):
    cols_df = pd.read_csv(cols_table, index_col=0)
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
        p1 = pairs_df.iloc[i, 0]
        p2 = pairs_df.iloc[i, 1]
        p1_str = p1.split('_')
        p2_str = p2.split('_')

        if p1_str[0] not in chr_df.index:
            # print p1_str[0]
            continue
        # if p1_str[0].find('chrM') != -1:
        #     print p1_str[0]

        if int(p1_str[1]) < int(p2_str[1]):
            distance = int(p2_str[1]) - int(p1_str[2])
        else:
            distance = int(p1_str[1]) - int(p2_str[2])
        if distance < 0:
            print p1, p2
            continue
        else:
            results[distance/100].add((p1, p2, distance))

    pairs = []
    pairs_set = set()

    for d in sorted(results.keys()):
        if len(results[d]) <= 10:
            for p in results[d]:
                pairs.append(p)
        else:
            candidates = list(results[d])
            candidates = random.sample(candidates, 10)
            for c in candidates:
                pairs.append(c)
                pairs_set.add(c[0])
                pairs_set.add(c[1])

    print "random pair selection complete!"

    f = open(signals_table, 'r')
    headers = f.readline()
    headers = headers.split(',')

    print len(headers)

    headers = headers[1:]

    print len(headers)

    headers_index = []
    for i in range(len(headers)):
        if headers[i] in pairs_set:
            headers_index.append(i)
            if headers[i].find('chrM') != -1:
                print headers[i]

    new_headers = [headers[i] for i in headers_index]
    new_pairs = set()
    for p in pairs:
        if p[0] in new_headers and p[1] in new_headers:
            new_pairs.add(p)

    pairs = new_pairs
    print len(pairs)

    print headers_index
    signals = []
    sample_id = 0
    for line in f.readlines():
        sample_id += 1
        line = line.split(',')
        line = line[1:]
        # print line[3813]
        # try:
            # for j in headers_index:
            #     try:
            #         n = float(line[j])
            #     except:
            #         print j, line[j], len(line[j]), headers[j]
        info = []
        for x in headers_index:
            try:
            # if line[x].isdigit():
                info.append(float(line[x]))
            except:
            # else:
                print line[x], headers[x], x, 'not digit?', sample_id
                info.append(0.0)
        signals.append(info)
        # except:
        #     print 'something is not number?!'

    signals_df = pd.DataFrame(signals)



    print 'chr9_131452500_131453270' in new_headers

    signals_df.columns = new_headers

    distance_corr_results = []
    # signals_df = pd.read_csv(signals_table, index_col=0)
    print "table load complete!"
    for p in pairs:
        r1 = signals_df[p[0]]
        r2 = signals_df[p[1]]
        cur_cor = np.corrcoef(r1, r2)[0, 1]
        # print r1[r1==0].shape[0]
        # print r1[pd.isnull(r1)].shape[0]
        # print r1
        print cur_cor, p[2]

        distance_corr_results.append((p[2], cur_cor))


    result_df = pd.DataFrame(distance_corr_results)
    #
    result_df.to_csv('peaks_distance_signal_correlation_100_10_2200.csv', index=None, header=None)
    return





close_peak_signals(['/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/' + x for x in
                    os.listdir('/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/')])

# close_peak_signals(['/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/' + x for x in os.listdir('/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/')])

# peaks_correlation('peaks_correlation.csv', 'peaks_correlation_columns.csv', None)
# print float(0)
# print float(format(0.00000, '.2f'))
# print float('')