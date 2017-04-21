"""
This module is to prepare the table the draw the cumulative curve
"""
import pandas as pd, numpy as np

def cumulative_cluster(stat_table, distance_table):
    """
    prepare the table for cumulative curve for each cluster
    :param stat_table:
    :param distance_table:
    :return:
    """
    stat_df = pd.read_csv(stat_table, sep='\t', index_col=0)
    distance_df = pd.read_csv(distance_table, index_col=0, header=None)
    distance_df.columns = ['TSS_distance']

    stat_df = stat_df.join(distance_df)

    clusters = stat_df['number of clusters'].unique().tolist()

    max_distance = stat_df['TSS_distance'].max()

    count_clusters = [0] * len(clusters)
    for k in clusters:
        count_clusters[k-1] = stat_df[stat_df['number of clusters'] == k].shape[0]
    # print count_clusters
    results = []
    for i in range(0, 100000, 500):
        cur_result = [0] * len(clusters)
        for j in clusters:
            cur_count = count_clusters[j-1]
            cur_df = stat_df[stat_df['number of clusters'] == j]
            percentage = cur_df[cur_df['TSS_distance']<=i].shape[0]*1.0/cur_count*100
            cur_result[j-1] = percentage
        # print cur_result
        results.append(cur_result)

    result_df = pd.DataFrame(results, columns=[x+1 for x in range(len(clusters))], index=np.asarray(range(0, 100000, 500))/1000)

    result_df.to_csv('clusters_distance.csv')

    results = []

    one_cluster = stat_df[(stat_df['number of clusters'] == 1)]
    # print pattern_df.shape
    more_than_one_df = stat_df[(stat_df['number of clusters'] != 1)]
    # print shape_df.shape

    for i in range(0, 100000, 500):
        one_percentage = one_cluster[one_cluster['TSS_distance'] <= i].shape[0] * 1.0 / one_cluster.shape[0] * 100
        not_one_percentage = more_than_one_df[more_than_one_df['TSS_distance'] <= i].shape[0] * 1.0 / more_than_one_df.shape[0] * 100
        results.append([one_percentage, not_one_percentage])

    result_df = pd.DataFrame(results, columns=['one', 'not_one'],
                             index=np.asarray(range(0, 100000, 500)) / 1000.0)

    result_df.to_csv('clusters_1_vs_not1_distance.csv')
    return result_df

def cumulative_pattern(stat_table, distance_table):
    """
    prepare the table for cumulative curve for pattern and shape change
    :param stat_table:
    :param distance_table:
    :return:
    """
    stat_df = pd.read_csv(stat_table, sep='\t', index_col=0)
    distance_df = pd.read_csv(distance_table, index_col=0, header=None)
    distance_df.columns = ['TSS_distance']

    stat_df = stat_df.join(distance_df)
    # print count_clusters
    results = []

    pattern_df = stat_df[(stat_df['BroadtoNarrow']==True) | (stat_df['Shift']==True) | \
                         (stat_df['ConcavetoConvex'] == True) | \
                         (stat_df['Pattern'] == True)]
    # print pattern_df.shape
    shape_df = stat_df[(stat_df['BroadtoNarrow']==False) & (stat_df['Shift']==False) & \
                       (stat_df['ConcavetoConvex'] == False) & \
                       (stat_df['Pattern'] == False)]
    # shape_df = stat_df[(stat_df['Other'] == True)]
    # print shape_df.shape

    for i in range(0, 100000, 500):
        pattern_percentage = pattern_df[pattern_df['TSS_distance'] <= i].shape[0] * 1.0 / pattern_df.shape[0] * 100
        shape_percentage = shape_df[shape_df['TSS_distance'] <= i].shape[0] * 1.0 / shape_df.shape[0] * 100
        results.append([pattern_percentage, shape_percentage])

    result_df = pd.DataFrame(results, columns=['Pattern', 'Shape'],
                             index=np.asarray(range(0, 100000, 500)) / 1000.0)

    result_df.to_csv('pattern_distance.csv')
    return result_df



# cumulative_cluster('./ref100/100_refregionstats.tsv', './ref100/total_region_distance.csv')

# cumulative_pattern('./ref100/100_refregionstats.tsv', './ref100/total_region_distance.csv')