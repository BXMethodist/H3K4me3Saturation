# This module is used to create a hashmap contains genome feature annotations with bin size 10.
# to facilitate the query of genome features for certain peaks.

import numpy as np, pandas as pd, pickle
from collections import defaultdict

def AnnotationToMap(annotation, outputmap):
    """
    :param annotation: gene annotation file downloaded for UCSC genome browser
    :param outputmap: the output map name
    :return: outputmap
    """
    df = pd.read_csv(annotation, sep='\t')

    print df.columns.values

    annotationmap = {}

    for i in range(df.shape[0]):
        gene_id = df.loc[i, '#name']
        gene_chr = df.loc[i, 'chrom']
        strand = 1 if df.loc[i, 'strand'] == '+' else -1
        tss = int(df.loc[i, 'txStart'])
        tts = int(df.loc[i, 'txEnd'])
        name2 = df.loc[i, 'name2']
        gene_obj = gene(gene_id, gene_chr, tss, tts, strand, name2)
        if gene_chr not in annotationmap:
            annotationmap[gene_chr] = defaultdict(set)

        annotationmap[gene_chr][tss/10*10].add(gene_obj)
        annotationmap[gene_chr][tts / 10 * 10].add(gene_obj)

    with open(outputmap + '.pkl', 'wb') as f:
        pickle.dump(annotationmap, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    return annotationmap

def feature_refmap(refmap, featuremap, cutoff1, cutoff2=0):
    refname = refmap[0:-4] if refmap.rfind('/') == -1 else refmap[refmap.rfind('/')+1:-4]
    featurename = featuremap[0:-4] if featuremap.rfind('/') == -1 else featuremap[featuremap.rfind('/')+1:-4]
    output_name = refname + featurename+str(cutoff1) + str(cutoff2)

    colnames = ['region_id', 'TSS', 'TSS_distance', 'gene_body', 'inter_gene', 'gene_name2']
    results = []
    best_results = []

    with open(refmap, 'rb') as f:
        refmap = pickle.load(f)
    f.close()

    with open(featuremap, 'rb') as f:
        featuremap = pickle.load(f)
    f.close()

    extend_range = max(cutoff1, cutoff2)

    for region in refmap:
        potential_genes = set()
        for key in range(region.start-extend_range, region.end+extend_range, 10):
            if key in featuremap[region.chromosome]:
                potential_genes = potential_genes.union(featuremap[region.chromosome][key])
        if len(potential_genes) == 0:
            results.append((region.id,
                            None,
                            None,
                            None,
                            True,
                            None))
            best_results.append((region.id, None, None, None, True, None))
            continue

        inter_gene = True

        best_distance = None
        best_overlap = 0
        best_gene_id = None
        best_type = None
        best_name2 = None

        for gene in potential_genes:
            if gene.strand == 1:
                tss_region_start, tss_region_end = gene.tss-cutoff1, gene.tss+cutoff2
                gene_body_region = (gene.tss+cutoff2, gene.tts) if gene.tss+cutoff2 < gene.tts else None
                if overlap(region.start, region.end, tss_region_start, tss_region_end) != 0:
                    results.append((region.id,
                                    gene.gene_id,
                                    (region.end-region.start)/2+region.start-gene.tss,
                                    None,
                                    None,
                                    gene.name2))
                    inter_gene = False

                    cur_overlap = overlap(region.start, region.end, tss_region_start, tss_region_end)

                    if cur_overlap > best_overlap:
                        best_distance = (region.end-region.start)/2+region.start-gene.tss
                        best_gene_id = gene.gene_id
                        best_overlap = cur_overlap
                        best_type = 'TSS'
                        best_name2 = gene.name2

                elif gene_body_region:
                    if overlap(region.start, region.end, gene_body_region[0], gene_body_region[1]):
                        results.append((region.id,
                                        None,
                                        None,
                                        gene.gene_id,
                                        None,
                                        gene.name2))
                        inter_gene = False

                        cur_overlap = overlap(region.start, region.end, gene_body_region[0], gene_body_region[1])

                        if cur_overlap > best_overlap:
                            best_distance = None
                            best_gene_id = gene.gene_id
                            best_overlap = cur_overlap
                            best_type = 'gene_body'
                            best_name2 = gene.name2

            elif gene.strand == -1:
                tss_region_start, tss_region_end = gene.tts - cutoff2, gene.tts + cutoff1
                gene_body_region = (gene.tss, gene.tts-cutoff2) if gene.tts - cutoff2 > gene.tss else None
                if overlap(region.start, region.end, tss_region_start, tss_region_end):
                    results.append((region.id,
                                    gene.gene_id,
                                    ((region.end - region.start) / 2 + region.start - gene.tts)*-1,
                                    None,
                                    None,
                                    gene.name2))
                    inter_gene = False

                    cur_overlap = overlap(region.start, region.end, tss_region_start, tss_region_end)

                    if cur_overlap > best_overlap:
                        best_distance = ((region.end - region.start) / 2 + region.start - gene.tts)*-1
                        best_gene_id = gene.gene_id
                        best_overlap = cur_overlap
                        best_type = 'TSS'
                        best_name2 = gene.name2

                elif gene_body_region:
                    if overlap(region.start, region.end, gene_body_region[0], gene_body_region[1]):
                        results.append((region.id,
                                        None,
                                        None,
                                        gene.gene_id,
                                        None,
                                        gene.name2))
                        inter_gene = False

                        cur_overlap = overlap(region.start, region.end, gene_body_region[0], gene_body_region[1])

                        if cur_overlap > best_overlap:
                            best_distance = None
                            best_gene_id = gene.gene_id
                            best_overlap = cur_overlap
                            best_type = 'gene_body'
                            best_name2 = gene.name2

        if inter_gene:
            results.append((region.id,
                            None,
                            None,
                            None,
                            True,
                            None))

        if best_type == 'TSS':
            best_results.append((region.id, best_gene_id, best_distance, None, None, best_name2))
        elif best_type == "gene_body":
            best_results.append((region.id, None, best_distance, best_gene_id, None, best_name2))
        else:
            best_results.append((region.id, None, None, None, True, None))

    df = pd.DataFrame(results, columns=colnames)
    df = df.set_index(['region_id'])
    df.to_csv(output_name+'.tsv', sep='\t')

    df2 = pd.DataFrame(best_results, columns=colnames)
    df2 = df2.set_index(['region_id'])
    df2.to_csv(output_name+'best_assign.tsv', sep='\t')

    return df

def combine_feature_cluster(feature, cluster):
    """
    :param feature: tsv with region feature
    :param cluster:  tsv with region cluster
    :return:
    """
    feature_df = pd.read_csv(feature, sep='\t', index_col=0)
    cluster_df = pd.read_csv(cluster, sep='\t', index_col=0)

    result_df = feature_df.join(cluster_df)

    result_df.to_csv(feature[:-4]+"with_cluster.tsv", sep='\t')

    return result_df

def StackedBarPlot(featuremap, groupby):
    """
    :param featuremap: dataframe
    :param groupby: feature groupby
    :return: percentage stacked bar plot
    """
    outputname = featuremap[:-4]+groupby+'stackedbarplot.tsv'

    df = pd.read_csv(featuremap, sep='\t', index_col=0)

    overall_counts = df.count(axis=0)

    results = [('All', overall_counts['TSS'], overall_counts['gene_body'], overall_counts['inter_gene'])]
    groups = df[groupby].unique()

    for group in groups:
        if not pd.isnull(group):
            sub_df = df[df[groupby]==group]
            sub_counts = sub_df.count(axis=0)
            results.append((group, sub_counts['TSS'], sub_counts['gene_body'], sub_counts['inter_gene']))
    featuremap_groups = pd.DataFrame(results, columns=['group', 'TSS', ' gene_body', 'inter_gene'])

    featuremap_groups.to_csv(outputname, sep='\t')

    return featuremap_groups

def overlap(start1, end1, start2, end2):
    """
    :param start1:
    :param end1:
    :param start2:
    :param end2:
    :return: return the overlap percentage for start1 and end1, if no overlap return 0
    """
    length = end1 - start1

    left = max(start1, start2)
    right = min(end1, end2)
    if left >= right:
        return 0
    else:
        return (right-left)*1.0/length

def comulative_TSS_plot(df, groupby, start=0, end=10000, bin=100):
    df = pd.read_csv(df, index_col=0, sep='\t')
    results = []
    colnames = ['distance', 'All'] + [group for group in df[groupby].unique() if not pd.isnull(group)]

    df_result = df[df['TSS'].notnull()]

    df_result['TSS_distance'] = df_result['TSS_distance'].abs()

    for i in range(start, end, bin):
        cur_result = []
        total = df_result.shape[0]
        sub_all = df_result[df_result['TSS_distance']>=start]
        sub_all = sub_all[sub_all['TSS_distance']<i]
        cur_result.append(i)
        cur_result.append(sub_all.shape[0]*1.0/total)
        for j in colnames[2:]:
            sub_set = sub_all[sub_all[groupby]==j]
            sub_set_total = df_result[df_result[groupby]==j]
            cur_result.append(sub_set.shape[0]*1.0/sub_set_total.shape[0])

        results.append(cur_result)

    return_df = pd.DataFrame(results, columns=colnames)
    return_df = return_df.set_index(['distance'])

    # return_df.plot()
    #
    # import matplotlib.pyplot as plt
    #
    # plt.show()


class gene():
    def __init__(self, gene_id, gene_chr, tss, tts, strand, name2):
        self.gene_id = gene_id
        self.gene_chr = gene_chr
        self.tss = tss
        self.tts = tts
        self.strand = strand
        self.name2 = name2





# feature_refmap('./pkl/75_combined_3kb.pkl', './pkl/hg38_RefSeq_allgene.pkl', 3000, 3000)
# combine_feature_cluster('75_combined_3kbhg38_RefSeq_allgene30003000.tsv', '75_combined_3kbstats.tsv')
#
# StackedBarPlot("75_combined_3kbhg38_RefSeq_allgene30003000with_cluster.tsv", "number of clusters")

# comulative_TSS_plot("75_combined_3kbhg38_RefSeq_allgene30003000with_cluster.tsv", "number of clusters")

