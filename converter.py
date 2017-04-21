"""
This module is used to convert reference map, gene annotation file to two kinds of structure of dictionary.
1. tier1 key: chr
tie2: Key, start and end.
value objects
2. tie1 chr
value: a sorted list of start, end.

then bisect could be used to find the closest features
"""
import os, pickle, pandas as pd
from collections import defaultdict
from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit


def AnnotationToMap(annotation, outputmap, tss_only):
    """ Create gene featurn annotation map, 1 lvl chromosome, 2 lvl location, 3 lvl gene obj
    :param annotation: gene annotation file downloaded for UCSC genome browser
    :param outputmap: the output map name
    :return: outputmap
    """
    df = pd.read_csv(annotation, sep='\t')

    print df.columns.values

    annotationmap = {}
    annotation_indexmap = {}

    if tss_only:
        outputmap += 'tss_only'

    for i in range(df.shape[0]):
        gene_id = df.loc[i, 'name']
        gene_chr = df.loc[i, 'chrom']
        strand = 1 if df.loc[i, 'strand'] == '+' else -1
        if strand == 1:
            tss = int(df.loc[i, 'txStart'])
            tts = int(df.loc[i, 'txEnd'])
        elif strand == -1:
            tts = int(df.loc[i, 'txStart'])
            tss = int(df.loc[i, 'txEnd'])
        name2 = df.loc[i, 'name2']
        gene_obj = gene(gene_id, gene_chr, tss, tts, strand, name2)
        if gene_chr not in annotationmap:
            annotationmap[gene_chr] = defaultdict(set)

        annotationmap[gene_chr][tss].add(gene_obj)
        if not tss_only:
            annotationmap[gene_chr][tts].add(gene_obj)

    for chromosome in annotationmap.keys():
        annotation_indexmap[chromosome] = sorted(annotationmap[chromosome].keys())

    with open(outputmap + '.pkl', 'wb') as f:
        pickle.dump(annotationmap, f, pickle.HIGHEST_PROTOCOL)
    f.close()

    with open(outputmap + '_index.pkl', 'wb') as f:
        pickle.dump(annotation_indexmap, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    return annotationmap, annotation_indexmap

def regiontomap(referencemap, outputmap):
    """
    Create reference map, 1 lvl chromosome, 2 lvl location, 3 lvl region obj
    :param annotation: reference map file
    :param outputmap: the output map name
    :return: outputmap
    """
    with open(referencemap, 'rb') as f:
        refmap = pickle.load(f)
    f.close()

    region_map = {}
    region_index_map = {}
    for region in refmap:
        if region.chromosome not in region_map:
            region_map[region.chromosome] = defaultdict(set)
        region_map[region.chromosome][region.start].add(region)
        region_map[region.chromosome][region.end].add(region)

    for chromosome in region_map.keys():
        region_index_map[chromosome] = sorted(region_map[chromosome].keys())
    
    with open(outputmap + '.pkl', 'wb') as f:
        pickle.dump(region_map, f, pickle.HIGHEST_PROTOCOL)
    f.close()

    with open(outputmap + '_index.pkl', 'wb') as f:
        pickle.dump(region_index_map, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    return region_map, region_index_map

def varianttomap(referencemap, outputmap):
    """
    Create reference map, 1 lvl chromosome, 2 lvl location, 3 lvl variant obj
    :param annotation: reference map file
    :param outputmap: the output map name
    :return: outputmap
    """
    with open(referencemap, 'rb') as f:
        refmap = pickle.load(f)
    f.close()

    variant_map = {}
    variant_index_map = {}
    for region in refmap:
        if region.chromosome not in variant_map:
            variant_map[region.chromosome] = defaultdict(set)
        for variant in region.variants:
            variant_map[region.chromosome][variant.left_boundary].add(variant)
            variant_map[region.chromosome][variant.right_boundary].add(variant)

    for chromosome in variant_map.keys():
        variant_index_map[chromosome] = sorted(variant_map[chromosome].keys())

    with open(outputmap + '.pkl', 'wb') as f:
        pickle.dump(variant_map, f, pickle.HIGHEST_PROTOCOL)
    f.close()

    with open(outputmap + '_index.pkl', 'wb') as f:
        pickle.dump(variant_index_map, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    return variant_map, variant_index_map

class gene():
    def __init__(self, gene_id, gene_chr, tss, tts, strand, name2):
        self.gene_id = gene_id
        self.gene_chr = gene_chr
        self.tss = tss
        self.tts = tts
        self.strand = strand
        self.name2 = name2


# AnnotationToMap('./pkl/hg19_RefSeq_refGene.txt', './pkl/hg19_RefSeq_refGene', True)



# './pkl/75_combined_3kb.pkl', './pkl/hg19_RefSeq_refGene.pkl'
# AnnotationToMap('./pkl/hg19_RefSeq_refGene.txt', './pkl/hg19_RefSeq_refGene', True)
# regiontomap('./ref100/100_refregion.pkl', './pkl/100_distance')
# varianttomap('./ref100/100_refregion.pkl', './pkl/100_variant_distance')
