import numpy as np
# from scipy.stats import pearsonr
# import pandas as pd
# from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
# from time import time
# import os

process = 8
chunks = []
chromosomes = [i for i in range(23)]
# print chromosomes

chunk_size = len(chromosomes) / process
reminder = len(chromosomes) % process

# print reminder

cur_index = 0

for i in range(process):
    if reminder > 0:
        chunks.append(chromosomes[cur_index + i * chunk_size:cur_index+(i + 1) * chunk_size + 1])
        cur_index += 1
        reminder -= 1
    else:
        chunks.append(chromosomes[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])


print chunks

# path = "/home/tmhbxx3/archive/test/"
#
# bed1 = "archive_tmhkxc48_BroadH3K4me3_test_ENCFF218WSN.chr21.Fnor.smooth.peaks.xls"
# bed2 = "archive_tmhkxc48_BroadH3K4me3_test_ENCFF621OIP.chr21.Fnor.smooth.peaks.xls"
#
#
# results = []
#
# for i in range(10, 310, 10):
#     cur_path = path+str(i)+"/pooled/"
#     df = pd.read_csv(cur_path+bed1, sep="\t")
#     # print df
#     peak_size = df.ix[:,"end"] - df.ix[:,"start"]
#     results.append((i, peak_size.values.T.mean()))
#
# result_file = open("./ENCFF218WSN.chr21_result_10to300_10.txt", "w")
# for x, y in results:
#     result_file.write(str(x) + "\t" + str(y) + "\n")
# result_file.close()
#
# results = []
#
# for i in range(10, 310, 10):
#     cur_path = path+str(i)+"/pooled/"
#     df = pd.read_csv(cur_path+bed2, sep="\t")
#     # print df
#     peak_size = df.ix[:,"end"] - df.ix[:,"start"]
#     results.append((i, peak_size.values.T.mean()))
#
# result_file = open("./ENCFF621OIP.chr21_result_10to300_10.txt", "w")
# for x, y in results:
#     result_file.write(str(x) + "\t" + str(y) + "\n")
# result_file.close()


# df = pd.read_csv("./csv/ENCFF621OIP.chr21_result_10to300_10.txt", sep="\t", index_col=0, header=None)
# df2 = pd.read_csv("./csv/ENCFF218WSN.chr21_result_10to300_10.txt", sep="\t", index_col=0, header=None)
# import matplotlib.pyplot as plt
#
#
#
# plot = df.plot()
#
# fig = plot.get_figure()
#
# fig.savefig('output.png')
#
#
# plot = df2.plot()
#
# fig = plot.get_figure()
#
# fig.savefig('output2.png')

#
# import matplotlib
#
# plt.figure()
# fig = df.plot()
#
# fig.savefig("./pictures/chr21", dpi=600, facecolor='w', edgecolor='w',
#                 orientation='portrait')
#
# peak_map = [0,1,0,0,0,1,1,0,0,1,0]
# peak_map = np.asarray(peak_map)
#
#
# sign = ((np.roll(peak_map, 1) - peak_map) != 0).astype(int)
# sign[0] = 0
# islandNumber = np.sum(sign)
# if peak_map[0] == 1:
#     islandNumber += 1
# if peak_map[-1] == 1:
#     islandNumber += 1
# print islandNumber/2

# import pandas as pd, numpy as np
# import matplotlib.pyplot as plt

# plt.xlim((0,340))
#
# df = pd.read_csv("chr3_clusters_count.csv", sep=",", index_col=0)
#
# result = {}
# for i in range(5):
#     result[i+1] ={}
#     for j in range(5,338):
#         result[i+1][j]=0
#
# for i in range(0, 5):
#     cur_values = df[df["#Cluster"]==i+1]['#sample_enriched'].values
#     for value in cur_values:
#         result[i+1][value]+=1
#     # print result[i+1]
#
#
#
# new_df = pd.DataFrame(result)
#
# new_df.columns =['1','2','3','4','5']
#
# for i in range(1,6):
#     new_df[str(i)] = new_df[str(i)].cumsum()/new_df[str(i)].sum()
#
#
#
# print new_df
#
# plot = new_df.plot(x=new_df.index, xlim=(0,340), colormap='prism', title='Number of Enriched Samples vs Number of Clusters')
# fig = plot.get_figure()
#
# fig.savefig("output.png")


# print df['#Sample_Clustered'].sum()


# df['cum_sum_sample_clustered'] = np.cumsum(df['#Sample_Clustered'].values)*1.0/df['#Sample_Clustered'].sum()
#
# import matplotlib.pyplot as plt
#
# plt.scatter(df['#Cluster'], df['#Sample_Clustered'], c=df['#Cluster'], cmap='prism', edgecolors='gray')
# plt.xlim((0,6))
# plt.ylim((0,340))
# plt.show()
