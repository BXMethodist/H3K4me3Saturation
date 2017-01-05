import numpy as np
from scipy.stats import pearsonr
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
from time import time
import os



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

peak_map = [0,1,0,0,0,1,1,0,0,1,0]
peak_map = np.asarray(peak_map)


sign = ((np.roll(peak_map, 1) - peak_map) != 0).astype(int)
sign[0] = 0
islandNumber = np.sum(sign)
if peak_map[0] == 1:
    islandNumber += 1
if peak_map[-1] == 1:
    islandNumber += 1
print islandNumber/2