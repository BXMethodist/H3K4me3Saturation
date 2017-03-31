import os, pandas as pd
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def distance(p0, p1, p2):
    """
    Calculate the largest distance bewteen a line and a point
    :param p0: first point to determine the line
    :param p1: second point to determine the line
    :param p2: this is the point that will be measured.
    :return: float number, distance
    """
    x0, y0 = p0
    x1, y1 = p1
    x2, y2 = p2
    nom = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1)
    denom = ((y2 - y1) ** 2 + (x2 - x1) ** 2) ** 0.5
    result = nom / denom
    return result

def cumulative_stat_vs_samplenumber_vs_cutoff(cutoffs, samplenumber_ranges, file_addresses, outputname):
    """
    Get the cumulative table for x: cutoffs, y: stat, different groups: sample_numbers
    :param cutoffs: peak height cutoffs
    :param file_addresses: file folder to keep the csv.
    :param outputname: output prefix.
    :return: 3 file names in the directory of outputname, region length file, region number file
    """
    coverages = pd.DataFrame(index=cutoffs, columns=[str(x) for x in samplenumber_ranges])
    length = pd.DataFrame(index=cutoffs, columns=[str(x) for x in samplenumber_ranges])
    region_number = pd.DataFrame(index=cutoffs, columns=[str(x) for x in samplenumber_ranges])

    if not file_addresses.endswith("/"):
        file_addresses += "/"
    for cutoff in cutoffs:
        df = pd.read_csv(file_addresses + str(cutoff) + '.csv', header=None)
        for r in samplenumber_ranges:
            coverages.ix[cutoff, str(r)] = df.iloc[r - 1, 1]
            length.ix[cutoff, str(r)] = df.iloc[r - 1, 2]
            region_number.ix[cutoff, str(r)] = df.iloc[r - 1, 3]

    coverages.to_csv(outputname+"cutoff_vs_coverage_samplenumber.csv")
    length.to_csv(outputname+"cutoff_vs_length_samplenumber.csv")
    region_number.to_csv(outputname+"cutoff_vs_regionnumber.csv")

    dfs = [coverages, length, region_number]
    outputnames = [outputname+"cutoff_vs_coverage_samplenumber.csv",
                   outputname + "cutoff_vs_length_samplenumber.csv",
                   outputname + "cutoff_vs_regionnumber.csv"]
    # for i in range(3):
    #     cur_df = dfs[i]
    #     ax = cur_df.plot(y=[str(x) for x in samplenumber_ranges], kind='line')
    #     cur_pic_name = outputnames[i]
    #     plt.tight_layout()
    #     plt.savefig(cur_pic_name[:-4]+".pdf", format='pdf', dpi=600, facecolor='w', edgecolor='w',
    #                 figsize=(3, 3))
    #     plt.close('all')
    return outputnames

def finalpoint_cutoff_vs_stat(cutoffs, file_addresses, number_sample, outputname, prefix):
    """
    generate table for x: cutoff, y: stat at final point
    :param cutoffs: peak height cutoffs
    :param file_addresses: file folder to keep the csv.
    :param number_sample: number of sample, this will determine the index use to get the final cumulative row.
    :param outputname: output prefix.
    :return: 1 file names in the directory of outputname
    """
    if not file_addresses.endswith("/"):
        file_addresses += "/"

    maps = [file_addresses + prefix + str(x) + ".csv" for x in cutoffs]

    table = []
    index = []

    # print maps
    for i in range(len(maps)):
        map = maps[i]
        df = pd.read_csv(map, header=None)
        info = df.iloc[number_sample-1, 1:]
        # print info
        table.append(info)
        # print map[:-4]
        index.append(cutoffs[i])

    df = pd.DataFrame(table, columns=None, index=index)
    # df.to_csv("q_vs_map.csv", header=None)

    df = df.sort_index()
    # print df

    # df = pd.read_csv("q_vs_map.csv", header=None, index_col=0)
    df.columns = ['Coverage (bp)', 'Length of Region', 'Number of Region (bp)']
    # df.to_csv(outputname + "cutoff_vs_map.csv")

    # for col in ['Coverage (bp)', 'Length of Region', 'Number of Region (bp)']:
    #     ax = df.plot(y=col, kind='line')
    #     plt.savefig(outputname + col+ "cutoff_vs_map" + ".pdf", format='pdf', dpi=600, facecolor='w', edgecolor='w',
    #                 figsize=(3, 3))
    #     plt.close('all')

    return df

def sample_number_vs_stat(cutoffs, file_addresses, sample_number, outputname):
    """
    get the table for plotting curves, x: number_sample, y: stat, group: cutoffs
    :param cutoffs: peak height cutoffs
    :param file_addresses: file folder to keep the csv.
    :param sample_number: used for generating number of rows
    :param outputname: output prefix.
    :return: 3 file names in the directory of outputname, region length file, region number file
    """
    if not file_addresses.endswith("/"):
        file_addresses += "/"

    cutoffs = sorted(cutoffs, reverse=True)

    coverage = pd.DataFrame(index=[i for i in range(1, sample_number+1)], columns=[str(i) for i in cutoffs])
    length = pd.DataFrame(index=[i for i in range(1, sample_number+1)], columns=[str(i) for i in cutoffs])
    region_number = pd.DataFrame(index=[i for i in range(1, sample_number+1)], columns=[str(i) for i in cutoffs])

    for cutoff in cutoffs:
        df = pd.read_csv(file_addresses+str(cutoff)+'.csv', header=None)
        coverage[str(cutoff)] = df.iloc[:, 1]
        length[str(cutoff)] = df.iloc[:, 2]
        region_number[str(cutoff)] = df.iloc[:, 3]

    coverage.to_csv(outputname+'samplenumber_coverage.csv')
    length.to_csv(outputname+'samplenumber_length.csv')
    region_number.to_csv(outputname+'samplenumber_regionnumber.csv')

    dfs = [coverage, length, region_number]
    outputnames = [outputname+'samplenumber_coverage.csv',
                   outputname + 'samplenumber_length.csv',
                   outputname + 'samplenumber_regionnumber.csv']
    # for i in range(3):
    #     cur_df = dfs[i]
    #     ax = cur_df.plot(y=[str(i) for i in cutoffs], kind='line')
    #     cur_pic_name = outputnames[i]
    #     plt.tight_layout()
    #     plt.savefig(cur_pic_name[:-4]+".pdf", format='pdf', dpi=600, facecolor='w', edgecolor='w',
    #                 figsize=(3, 3))
    #     plt.close('all')
    return outputnames