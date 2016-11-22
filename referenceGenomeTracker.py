import pandas as pd, numpy as np
import os



genome = {}

file = open("hg19_chr_sizes.txt", "r")
for line in file.readlines():
    chrName, chrSize = line.split("\t")
    genome[chrName] = np.zeros(int(chrSize.rstrip()))

file.close()

df = pd.read_excel("regions.xls")
# print df.ix[0,:]['chr']

rowNum = df.shape[0]

for i in range(len(rowNum)):
    start = df.ix[i, :]['start']
    end = df.ix[i, :]['end']
    chrName = df.ix[i, :]['chr']


