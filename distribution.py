"""
This module is used to calculate the distribution
"""


import pandas as pd
def get_distribution(file_path, step=100):
    """
    this is function take a txt file, with single columns and return is distribution based on step
    :param file_path:
    :param step:
    :return:
    """
    df = pd.read_csv(file_path, index_col=None, header=None)
    df.columns = ['distance']
    indexes = df['distance'].unique().tolist()
    indexes = sorted(indexes, key=lambda x: int(x))
    print len(indexes)
    results = []
    counts = {}
    cum = 0
    for i in indexes:
        cur_df = df[df['distance']==i]
        counts[i] = cur_df.shape[0]
        cum += counts[i]
        p = cum*1.0/df.shape[0]
        results.append([i, p])
        # if p > 0.99:
        #     break
        # print 'count is done'

    # max_distance = df['distance'].max()
    #
    # step = 10
    # n = 1
    # m = 1
    # k = 1
    # while n * step < max_distance:
    #     cur_len = df[df.ix[:, 0] < n * step].shape[0]
    #     results.append([n*step, cur_len * 1.0 / df.shape[0]])
    #     if cur_len*1.0/df.shape[0] > 0.99:
    #         break
    #     n += 1
    #     m += 1
    #     if m >= 100:
    #         step +=10
    #         m = 1

    # for i in indexes:
    #     cur_len = df[df.ix[:, 0] < i].shape[0]
    #     results.append([i, cur_len*1.0/df.shape[0]])
    #     if cur_len*1.0/df.shape[0] > 0.99:
    #         break
    result_df = pd.DataFrame(results)
    result_df.to_csv(file_path[:-4]+'distribution.csv', index=None, header=None)
    return result_df

# get_distribution('./extend_100_10/gene_up_distance_distribution.txt')
# get_distribution('./extend_100_10/gene_down_distance_distribution.txt')
# get_distribution('./extend_100_10/gene_total_distance_distribution.txt')
# get_distribution('./extend_100_10/region_up_distance_distribution.txt')
# get_distribution('./extend_100_10/region_down_distance_distribution.txt')
# get_distribution('./extend_100_10/region_total_distance_distribution.txt')
# get_distribution('./ref100/regions_per_samples.txt', step=1)
# get_distribution('./extend_100_10_2200/100_10_2200_extend_merge_distance.txt', step=100)

# get_distribution('/Users/boxia/Desktop/reference_map/QC/bowtie_results_ENC_percentage.txt', step=0.01)
# get_distribution('/Users/boxia/PycharmProjects/H3K4me3Saturation/extend_100_10_2200/peaks_distance_signal_correlation_100.csv')
# get_distribution('extend_distance_percentage.txt', step=10)

# get_distribution('region_gap_distance.txt', step=100)

# get_distribution('100_refmap_peaks_distance.txt', step=100)

# get_distribution('gene_to_gene_tss_distance.txt', step=100)

# df = pd.read_csv('100_refmap_peaks_distance.txt', index_col=None, header=None)


#
# df.columns=['distance']
#
# distance = df['distance'].values
#
# distance = distance[distance<1000]
#
# import matplotlib.pyplot as plt
#
# plt.hist(distance, bins=100)
# plt.show()