# This module is used to create a hashmap contains genome feature annotations with bin size 10.
# to facilitate the query of genome features for certain peaks.

# 1. Convert and Annotation file to pkl by AnnotationToMap
# 2. Map annotation to peaks by 'feature_refmap'
# 3. combine featured referencemap with cluster_category map by 'combine_feature_cluster'
# 4. Generate excel table for plot by 'StackedBarPlot'

import numpy as np, pandas as pd, pickle
from collections import defaultdict
from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit

def AnnotationToMap(annotation, outputmap):
    """ Create gene featurn annotation map, 1 lvl chromosome, 2 lvl location, 3 lvl gene obj
    :param annotation: gene annotation file downloaded for UCSC genome browser
    :param outputmap: the output map name
    :return: outputmap
    """
    df = pd.read_csv(annotation, sep='\t')

    print df.columns.values

    annotationmap = {}

    for i in range(df.shape[0]):
        gene_id = df.loc[i, 'name']
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

def feature_refmap(refmap, featuremap, cutoff1, cutoff2=0, outputdir='./pkl/'):
    """
    Map the refmap with feature map based on the cutoff, cutoff1 is for -, cutoff2 is for +,
    for example, -3000, 3000 means map feature by +- 3000 bp
    :param refmap:
    :param featuremap:
    :param cutoff1:
    :param cutoff2:
    :param outputdir:
    :return:
    """
    refname = refmap[0:-4] if refmap.rfind('/') == -1 else refmap[refmap.rfind('/')+1:-4]
    featurename = featuremap[0:-4] if featuremap.rfind('/') == -1 else featuremap[featuremap.rfind('/')+1:-4]
    output_name = refname + featurename+str(cutoff1) + str(cutoff2)

    colnames = ['variant_id', 'region_id', 'TSS', 'TSS_distance', 'gene_body', 'inter_gene', 'gene_name2']
    results = []
    best_results = []

    with open(refmap, 'rb') as f:
        refmap = pickle.load(f)
    f.close()

    with open(featuremap, 'rb') as f:
        featuremap = pickle.load(f)
    f.close()

    extend_range = max(cutoff1, cutoff2)

    count_gene_promoter = set()
    count_gene_body = set()
    count_transcript_promoter = set()
    count_transcript_body = set()


    for region in refmap:
        for variant in region.variants:
            potential_genes = set()
            left = variant.start
            right = variant.end
            for key in range(left-extend_range, right+extend_range, 10):
                if key in featuremap[region.chromosome]:
                    potential_genes = potential_genes.union(featuremap[region.chromosome][key])
            if len(potential_genes) == 0:
                results.append((variant.id,
                                region.id,
                                None,
                                None,
                                None,
                                True,
                                None))
                best_results.append((variant.id, None, None, None, True, None))
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
                    if overlap(variant.start, variant.end, tss_region_start, tss_region_end) != 0:
                        results.append((variant.id,
                                        region.id,
                                        gene.gene_id,
                                        variant.end-variant.start-gene.tss,
                                        None,
                                        None,
                                        gene.name2))
                        inter_gene = False

                        count_transcript_promoter.add(gene.gene_id)
                        count_gene_promoter.add(gene.name2)
                        if gene.gene_id in count_transcript_body:
                            count_transcript_body.remove(gene.gene_id)
                        if gene.name2 in count_gene_body:
                            count_gene_body.remove(gene.name2)

                        cur_overlap = overlap(left, right, tss_region_start, tss_region_end)

                        if cur_overlap > best_overlap:
                            best_distance = variant.end-variant.start-gene.tss
                            best_gene_id = gene.gene_id
                            best_overlap = cur_overlap
                            best_type = 'TSS'
                            best_name2 = gene.name2

                    elif gene_body_region:
                        if overlap(left, right, gene_body_region[0], gene_body_region[1]):
                            results.append((variant.id,
                                            region.id,
                                            None,
                                            None,
                                            gene.gene_id,
                                            None,
                                            gene.name2))
                            inter_gene = False

                            count_transcript_body.add(gene.gene_id)
                            count_gene_body.add(gene.name2)

                            cur_overlap = overlap(left, right, gene_body_region[0], gene_body_region[1])

                            if cur_overlap > best_overlap:
                                best_distance = None
                                best_gene_id = gene.gene_id
                                best_overlap = cur_overlap
                                best_type = 'gene_body'
                                best_name2 = gene.name2

                elif gene.strand == -1:
                    tss_region_start, tss_region_end = gene.tts - cutoff2, gene.tts + cutoff1
                    gene_body_region = (gene.tss, gene.tts-cutoff2) if gene.tts - cutoff2 > gene.tss else None
                    if overlap(left, right, tss_region_start, tss_region_end):
                        results.append((variant.id,
                                        region.id,
                                        gene.gene_id,
                                        (variant.end-variant.start - gene.tts)*-1,
                                        None,
                                        None,
                                        gene.name2))
                        inter_gene = False

                        count_transcript_promoter.add(gene.gene_id)
                        count_gene_promoter.add(gene.name2)
                        if gene.gene_id in count_transcript_body:
                            count_transcript_body.remove(gene.gene_id)
                        if gene.name2 in count_gene_body:
                            count_gene_body.remove(gene.name2)

                        cur_overlap = overlap(left, right, tss_region_start, tss_region_end)

                        if cur_overlap > best_overlap:
                            best_distance = (variant.end-variant.start - gene.tts)*-1
                            best_gene_id = gene.gene_id
                            best_overlap = cur_overlap
                            best_type = 'TSS'
                            best_name2 = gene.name2

                    elif gene_body_region:
                        if overlap(left, right, gene_body_region[0], gene_body_region[1]):
                            results.append((variant.id,
                                            region.id,
                                            None,
                                            None,
                                            gene.gene_id,
                                            None,
                                            gene.name2))
                            inter_gene = False

                            count_transcript_body.add(gene.gene_id)
                            count_gene_body.add(gene.name2)

                            cur_overlap = overlap(left, right, gene_body_region[0], gene_body_region[1])

                            if cur_overlap > best_overlap:
                                best_distance = None
                                best_gene_id = gene.gene_id
                                best_overlap = cur_overlap
                                best_type = 'gene_body'
                                best_name2 = gene.name2

            if inter_gene:
                results.append((variant.id,
                                region.id,
                                None,
                                None,
                                None,
                                True,
                                None))

            if best_type == 'TSS':
                best_results.append((variant.id, region.id, best_gene_id, best_distance, None, None, best_name2))
            elif best_type == "gene_body":
                best_results.append((variant.id, region.id, None, best_distance, best_gene_id, None, best_name2))
            else:
                best_results.append((variant.id, region.id, None, None, None, True, None))

    outputdir = outputdir +'/' if not outputdir.endswith('/') else outputdir

    df = pd.DataFrame(results, columns=colnames)
    df = df.set_index(['variant_id'])
    df.to_csv(outputdir+output_name+'.tsv', sep='\t')

    df2 = pd.DataFrame(best_results, columns=colnames)
    df2 = df2.set_index(['variant_id'])
    df2.to_csv(outputdir+output_name+'best_assign.tsv', sep='\t')

    print "there is ", len(count_gene_promoter), " gene promoter covered by the reference map"
    print "there is ", len(count_gene_body), " gene body covered by the reference map and promoter not covered by " \
                                             "reference map"
    print "there is ", len(count_transcript_promoter), " transcript promoter covered by the reference map"
    print "there is ", len(count_transcript_body), " transcript body covered by the reference map and promoter not covered by " \
                                             "reference map"


    return df

def combine_feature_cluster(feature, cluster):
    """
    Combine the features for region (TSS, intergene...) with number of clusters
    :param feature: tsv with region feature
    :param cluster:  tsv with region cluster
    :return: combine gene feature (TSS, genebody...) information with cluster, category(number of clusters, BN, PT..)
    information, by using variant id as index
    """
    feature_df = pd.read_csv(feature, sep='\t', index_col=0)
    cluster_df = pd.read_csv(cluster, sep='\t', index_col=0)

    result_df = feature_df.join(cluster_df)

    result_df.to_csv(feature[:-4]+"with_cluster.tsv", sep='\t')

    return result_df

def StackedBarPlot(featuremap, groupby, groupby2=[]):
    """
    :param featuremap: dataframe
    :param groupby: feature groupby, if only one features in groupby, it will group the different value in that
    features towards TSS, genebody....
    if more than one feature in groupby, it will group different feature with TSS, gene body.
    if groupby feature is 1 and groupby2 feature is not empty, it means get the number of feature group 2 groupby feature group 1.
    :return: a dataframe used for draw the percentage stacked bar plot
    """
    if len(groupby) == 1 and len(groupby2) == 0:
        groupby = groupby[0]
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
    elif len(groupby) > 1:
        outputname = featuremap[:-4] + '_'.join(groupby) + 'stackedbarplot.tsv'
        df = pd.read_csv(featuremap, sep='\t', index_col=0)
        overall_counts = df.count(axis=0)
        results = [('All', overall_counts['TSS'], overall_counts['gene_body'], overall_counts['inter_gene'])]

        for group in groupby:
            sub_df = df[df[group]==True]
            sub_counts = sub_df.count(axis=0)
            results.append((group, sub_counts['TSS'], sub_counts['gene_body'], sub_counts['inter_gene']))
        featuremap_groups = pd.DataFrame(results, columns=['group', 'TSS', ' gene_body', 'inter_gene'])

        featuremap_groups.to_csv(outputname, sep='\t')

    elif len(groupby) == 1 and len(groupby2) != 0:
        outputname = featuremap[:-4] + '_'.join(groupby) + "_".join(groupby2) + 'stackedbarplot.tsv'
        df = pd.read_csv(featuremap, sep='\t', index_col=0)
        results = []

        for feature in groupby2:
            df.loc[df[feature]==False, feature] = np.nan

        for group in df[groupby[0]].unique():
            if not pd.isnull(group):
                sub_df = df[df[groupby[0]]==group]
                sub_counts = sub_df.count(axis=0)
                cur_result =[group] + [sub_counts[feature] for feature in groupby2]
                results.append(cur_result)
        featuremap_groups = pd.DataFrame(results, columns=['group']+ [feature for feature in groupby2])

        featuremap_groups.to_csv(outputname, sep='\t')




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
    """
    This function is used to check the distance of each cluster towards TSS
    :param df:
    :param groupby:
    :param start:
    :param end:
    :param bin:
    :return:
    """
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

    return_df.plot()

    import matplotlib.pyplot as plt

    plt.show()


class gene():
    def __init__(self, gene_id, gene_chr, tss, tts, strand, name2):
        self.gene_id = gene_id
        self.gene_chr = gene_chr
        self.tss = tss
        self.tts = tts
        self.strand = strand
        self.name2 = name2



# AnnotationToMap('./pkl/hg19_RefSeq_refGene.txt', './pkl/hg19_RefSeq_refGene')

feature_refmap('./pkl/75_combined_3kb.pkl', './pkl/hg19_RefSeq_refGene.pkl', 3000, 3000, outputdir='./pkl')
# combine_feature_cluster('./pkl/75_combined_3kbhg19_RefSeq_refGene30003000.tsv', './pkl/75_combined_3kbstats.tsv')
# # #
#StackedBarPlot("./pkl/75_combined_3kbhg19_RefSeq_refGene30003000with_cluster.tsv",['Pattern'])
               #['number of clusters'])
               #['BroadtoNarrow',	'ConcavetoConvex',	'Shift',	'Pattern',	'Other'])

# comulative_TSS_plot("./pkl/75_combined_3kbhg19_RefSeq_refGene30003000with_cluster.tsv", "number of clusters")

