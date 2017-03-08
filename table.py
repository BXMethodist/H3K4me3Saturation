import os, pandas as pd
import numpy as np


def distance(p0, p1, p2):  # p3 is the point
    x0, y0 = p0
    x1, y1 = p1
    x2, y2 = p2
    nom = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1)
    denom = ((y2 - y1) ** 2 + (x2 - x1) ** 2) ** 0.5
    result = nom / denom
    return result


maps = [x for x in os.listdir("./") if x.endswith(".csv") and x.find("refmap") == -1 and x != 'q_vs_map.csv']

table = []
index = []

for map in maps:
    if int(map[:-4]) != 25 and int(map[:-4]) >= 10:
        df = pd.read_csv(map, header=None)
        info = df.iloc[335, 1:]
        # print info
        table.append(info)
        # print map[:-4]
        index.append(int(map[:-4]))

df = pd.DataFrame(table, columns=None, index=index)
# df.to_csv("q_vs_map.csv", header=None)

df = df.sort_index()
# print df

# df = pd.read_csv("q_vs_map.csv", header=None, index_col=0)
df.columns = ['Coverage', 'Length of Region', 'Number of Region']

df['Number of Region'] = np.log(df['Number of Region'])
df['Coverage'] = np.log(df['Coverage'])

max_dis = 0
cutoff = None
for i in range(10, 310, 10):
    cur_dis = distance((10, df.iloc[0, 0]),
                       (300, df.iloc[-1, 0]),
                       (i * 10, df.iloc[i, 0]))
    if cur_dis > max_dis:
        max_dis = cur_dis
        cutoff = None
print cutoff

max_dis = 0
cutoff = None
for i in range(10, 310, 10):
    cur_dis = distance((10, df.iloc[0, 2]),
                       (300, df.iloc[-1, 2]),
                       (i * 10, df.iloc[i, 2]))
    if cur_dis > max_dis:
        max_dis = cur_dis
        cutoff = None
print cutoff

col = 'Number of Region'
ax = df.plot(y=col)
ax.set_xlim((10, 300))
ax.set_xlabel('Peak Cutoff')
ax.set_ylabel(col)

import matplotlib.pyplot as plt

plt.show()



