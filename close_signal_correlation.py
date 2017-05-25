import os, pandas as pd, numpy as np, pickle, random, bisect
from refMapUtils import load_obj

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
            results.append(float(format(np.sum(cur_signals), '.2f')))
        else:
            results.append(None)
    return results

def close_peak_signals(peak_tables):
    sample_peak_numbers = {}
    total_number = 0
    sample_peaks = {}
    sample_peaks_list = {}

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

    random_size = 10000

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

    for wig in wigs:
        wigs_signals.append(get_wig_signals(wig, correlations_regions))


    df = pd.DataFrame(wigs_signals)
    df.columns= ['_'.join(list([str(z) for z in region])) for region in correlations_regions]

    df.to_csv('peaks_correlation.csv')

    df_cols = pd.DataFrame(correlations_pairs)

    df_cols.to_csv('peaks_correlation_columns.csv')

    return df, df_cols

def peaks_correlation(signals_table, cols_table, random_pairs):
    cols_df = pd.read_csv(cols_table, index_col=0)
    if random_pairs is not None:
        pairs_df = cols_df.sample(random_pairs)
    else:
        pairs_df = cols_df.copy()
    results = []
    for i in range(pairs_df.shape[0]):
        p1 = pairs_df.iloc[i, 0]
        p2 = pairs_df.iloc[i, 1]
        p1_str = p1.split('_')
        p2_str = p2.split('_')

        if int(p1_str[1]) < int(p2_str[1]):
            distance = int(p2_str[1]) - int(p1_str[2])
        else:
            distance = int(p1_str[1]) - int(p2_str[2])
        if distance < 0:
            print p1, p2
            continue
        else:
            results.append((p1, p2, distance))

    print "random pair selection complete!"

    distance_corr_results = []
    signals_df = pd.read_csv(signals_table, index_col=0)
    for result in results:
        r1 = signals_df[result[0]]
        r2 = signals_df[result[1]]
        cur_cor = np.corrcoef(r1, r2)[0,1]
        distance_corr_results.append((result[2], cur_cor))

    result_df = pd.DataFrame(distance_corr_results)

    result_df.to_csv('peaks_distance_signal_correlation_100.csv', index=None, header=None)
    return





# close_peak_signals(['/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/' + x for x in
#                     os.listdir('/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/')])

# close_peak_signals(['/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/' + x for x in os.listdir('/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/')])

peaks_correlation('peaks_correlation.csv', 'peaks_correlation_columns.csv', None)
