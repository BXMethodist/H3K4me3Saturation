"""
This module is used to collect genes list from analysis results and run the david path way analysis
"""

import os, pandas as pd, numpy as np

def DavidGenes(table, features, logic_operater=None, sep='\t', index_col=0):
    """
    get the genelist by featuren and featuren_value
    :param table: table file path
    :param feature: a dictionary contain feature name, (column name) and feature value as pair.
    :param logic_operater: None, or, and
    :return:
    """
    df = pd.read_csv(table, sep=sep, index_col=index_col)
    if logic_operater is None:
        for feature, feature_value in features.items():
            sub_df = df[df[feature] == feature_value]
            gene_list = sub_df['gene_name2'].unique()
            result_file = open(feature+"_gene_list.txt", "w")
            print gene_list
            for gene in gene_list:
                result_file.write(str(gene)+'\n')
            result_file.close()
    if logic_operater == 'or':
        selections = np.asarray([False] * df.shape[0])
        for feature, feature_value in features.items():
            selections = selections | df[feature] == feature_value
        sub_df = df[selections]
        gene_list = sub_df['gene_name2'].unique()
        combined_name = "_".join(features.keys())
        result_file = open(combined_name+"_gene_list.txt", "w")
        for gene in gene_list:
            result_file.write(str(gene)+'\n')
        result_file.close()
    if logic_operater == 'and':
        selections = np.asarray([True] * df.shape[0])
        for feature, feature_value in features.items():
            selections = selections & df[feature] == feature_value
        sub_df = df[selections]
        gene_list = sub_df['gene_name2'].unique()
        combined_name = "_".join(features.keys())
        result_file = open(combined_name+"_gene_list.txt", "w")
        for gene in gene_list:
            result_file.write(str(gene)+'\n')
        result_file.close()
    return

def RandomGenes(annotation_file, gene_name_column, gene_number):
    """
    :param annotation_file: annotation file
    :param gene_name_column: column_name for official gene symbol
    :param gene_number: number of genes in random list
    :return:
    """
    df = pd.read_csv(annotation_file, sep='\t')
    gene_list = df['name2'].unique()
    random_list = np.random.choice(gene_list, replace=False, size=gene_number)
    result_file = open("random_gene_list.txt", "w")
    for gene in random_list:
        result_file.write(str(gene) + '\n')
    result_file.close()
    return result_file

def KEGG_to_csv(david_lists):
    """ convert david KEGG table to csv, remove non-necessary columns
    :param david_lists:
    :return: a list of tuples (df, df_name)
    """
    columns =["Term", "Count", "Genes", "Fold Enrichment", "Benjamini"]
    dfs = []
    for l in david_lists:
        df = pd.read_csv(l, sep='\t')
        new_terms = [term[1] for term in df['Term'].str.split(":")]
        df["Term"] = new_terms
        df = df[columns]
        df = df.set_index(["Term"])
        df.to_csv(l[:-4]+"_converted.csv")
        df_name = l[l.rfind("/")+1:-4]
        dfs.append((df, df_name))
    result_df = MergeDavid(dfs)
    result_df.to_csv("merged_KEGG.csv")
    return result_df

def MergeDavid(david_df_lists):
    """
    Used to merge different David lists
    :param david_df_lists:a list of tuples (df, df_name)
    :return:
    """
    result_df = david_df_lists[0][0]
    for df, df_name in david_df_lists[1:]:
        result_df = result_df.join(df, rsuffix="_"+df_name)
    return result_df





# file_path = './pkl/75_combined_3kbhg19_RefSeq_refGene30003000with_cluster.tsv'
# #
# DavidGenes(file_path, {"BroadtoNarrow":True, "Shift":True, "Pattern":True, "ConcavetoConvex":True}, logic_operater="or")

# "BroadtoNarrow":True, "Shift":True, "Pattern":True, "ConcavetoConvex":True,

dfs = ['./pkl/KEGG_SH_BN_CC_PT.txt', './pkl/KEGG_Other_all.txt','./pkl/KEGG_Random_all.txt',]
KEGG_to_csv(dfs)