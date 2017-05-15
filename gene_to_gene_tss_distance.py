import os, pandas as pd
from itertools import combinations
import numpy as np
from collections import defaultdict

def tss_distance(comb, df):
    name1, name2 = comb
    tss1 = []
    tss2 = []

    df1 = df[df['name2'] == name1]

    for i in range(df1.shape[0]):
        if df1.iloc[i]['strand'] == '+':
            tss1.append(int(df1.iloc[i]['txStart']))
        elif df1.iloc[i]['strand'] == '-':
            tss1.append(int(df1.iloc[i]['txEnd']))

    df2 = df[df['name2'] == name2]

    for i in range(df2.shape[0]):
        if df2.iloc[i]['strand'] == '+':
            tss2.append(int(df2.iloc[i]['txStart']))
        elif df2.iloc[i]['strand'] == '-':
            tss2.append(int(df2.iloc[i]['txEnd']))

    best_distance = None
    for i in range(len(tss1)):
        for j in range(len(tss2)):
            if best_distance is None or best_distance > abs(tss1[i] - tss2[j]):
                best_distance = abs(tss1[i] - tss2[j])
    return best_distance


gtf_df = pd.read_csv('./pkl/hg19_RefSeq_refGene.txt', sep='\t', index_col=0)

# print gtf_df.columns

distances = []

for chrom in gtf_df['chrom'].unique():
    tsses = defaultdict(set)
    names = defaultdict(set)

    cur_df = gtf_df[gtf_df['chrom']==chrom].copy()

    for i in range(cur_df.shape[0]):
        if cur_df.iloc[i]['strand'] == '+':
            tsses[int(cur_df.iloc[i]['txStart'])].add(cur_df.iloc[i]['name2'])
            names[cur_df.iloc[i]['name2']].add(int(cur_df.iloc[i]['txStart']))
        elif cur_df.iloc[i]['strand'] == '-':
            tsses[int(cur_df.iloc[i]['txEnd'])].add(cur_df.iloc[i]['name2'])
            names[cur_df.iloc[i]['name2']].add(int(cur_df.iloc[i]['txEnd']))





    for gene in cur_df['name2'].unique():
        df1 = cur_df[cur_df['name2'] == gene]
        tss = []
        for i in range(df1.shape[0]):
            if df1.iloc[i]['strand'] == '+':
                tss.append(int(df1.iloc[i]['txStart']))
            elif df1.iloc[i]['strand'] == '-':
                tss.append(int(df1.iloc[i]['txEnd']))
        tsses.append(np.mean(tss))
    tsses = sorted(tsses)
    for i in range(len(tsses) -1):
        distances.append(tsses[i+1]-tsses[i])

    # genes = combinations(cur_df['name2'].unique(), 2)
    #
    # for comb in genes:
    #     distances.append(tss_distance(comb, gtf_df))
#         break
#
# print distances
f = open('gene_to_gene_tss_distance.txt', 'w')

for line in distances:
    f.write(str(line)+'\n')
f.close()

# for name2 in gtf_df['name2'].unique():
