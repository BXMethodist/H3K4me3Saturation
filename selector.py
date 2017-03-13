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

    colnames = ['region_id', 'TSS', 'TSS_distance', 'gene_body', 'inter_gene']
    results = []

    with open(refmap, 'rb') as f:
        refmap = pickle.load(f)
    f.close()

    with open(featuremap, 'rb') as f:
        featuremap = pickle.load(f)
    f.close()

    for region in refmap:
        potential_genes = set()
        for key in range(region.start-cutoff, region.end+cutoff, 10):
            if key in featuremap[region.chromosome]:
                potential_genes = potential_genes.union(featuremap[region.chromosome][key])
        if len(potential_genes) == 0:
            results.append((region.id,
                            None,
                            None,
                            None,
                            True))

        for gene in potential_genes:
            inter_gene = True
            if gene.strand == 1:
                tss_region = (gene.tss-cutoff, gene.tss+cutoff)
                gene_body_region = (gene.tss+cutoff, gene.tts) if gene.tss+cutoff > gene.tss else None
                if (region.start <= tss_region[0] and region.end >= tss_region[1]) or \
                        (tss_region[0]<=region.start and tss_region[1]>=region.end) or \
                        (tss_region[0]<=region.start and tss_region[1]>=region.start) or \
                        (tss_region[0]<=region.end and tss_region[1]>=region.end):
                    results.append((region.id,
                                    gene.gene_id,
                                    (region.end-region.start)/2+region.start-gene.tss,
                                    None,
                                    None))
                    inter_gene = False
                elif gene_body_region:
                    if (region.start <= gene_body_region[0] and region.end >= gene_body_region[1]) or \
                            (gene_body_region[0] <= region.start and gene_body_region[1] >= region.end) or \
                            (gene_body_region[0] <= region.start and gene_body_region[1] >= region.start) or \
                            (gene_body_region[0] <= region.end and gene_body_region[1] >= region.end):
                        results.append((region.id,
                                        None,
                                        None,
                                        gene.gene_id,
                                        None))
                        inter_gene = False
                if inter_gene:
                    results.append((region.id,
                                    None,
                                    None,
                                    None,
                                    True))
            elif gene.strand == -1:
                tss_region = (gene.tts - cutoff, gene.tts + cutoff)
                gene_body_region = (gene.tss, gene.tts-cutoff) if gene.tss + cutoff > gene.tss else None
                if (region.start <= tss_region[0] and region.end >= tss_region[1]) or \
                        (tss_region[0] <= region.start and tss_region[1] >= region.end) or \
                        (tss_region[0] <= region.start and tss_region[1] >= region.start) or \
                        (tss_region[0] <= region.end and tss_region[1] >= region.end):
                    results.append((region.id,
                                    gene.gene_id,
                                    ((region.end - region.start) / 2 + region.start - gene.tts)*-1,
                                    None,
                                    None))
                    inter_gene = False
                elif gene_body_region:
                    if (region.start <= gene_body_region[0] and region.end >= gene_body_region[1]) or \
                            (gene_body_region[0] <= region.start and gene_body_region[1] >= region.end) or \
                            (gene_body_region[0] <= region.start and gene_body_region[1] >= region.start) or \
                            (gene_body_region[0] <= region.end and gene_body_region[1] >= region.end):
                        results.append((region.id,
                                        None,
                                        None,
                                        gene.gene_id,
                                        None))
                        inter_gene = False
                if inter_gene:
                    results.append((region.id,
                                    None,
                                    None,
                                    None,
                                    True))

    df = pd.DataFrame(results, columns=colnames)
    df = df.set_index(['region_id'])
    df.to_csv(output_name+'.tsv', sep='\t')
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


class gene():
    def __init__(self, gene_id, gene_chr, tss, tts, strand, name2):
        self.gene_id = gene_id
        self.gene_chr = gene_chr
        self.tss = tss
        self.tts = tts
        self.strand = strand
        self.name2 = name2


# feature_refmap('./pkl/75_combined_3kb.pkl', './pkl/hg38_RefSeq_allgene.pkl', 3000)
# combine_feature_cluster('75_combined_3kbhg38_RefSeq_allgene3000.tsv', '75_combined_3kbstats.tsv')

StackedBarPlot("75_combined_3kbhg38_RefSeq_allgene3000with_cluster.tsv", "number of clusters")