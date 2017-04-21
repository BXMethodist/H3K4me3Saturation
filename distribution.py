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
    results = []
    biggest = df.max().values[0]
    print biggest
    n = 1
    while n * step < biggest:
        cur_len = df[df.ix[:, 0] < n*step].shape[0]
        n += 1
        results.append([cur_len*1.0/df.shape[0]])
    result_df = pd.DataFrame(results)
    result_df.to_csv(file_path[:-4]+'distribution.csv', index=None, header=None)
    return result_df

# get_distribution('./original_raw_peak.csv')

