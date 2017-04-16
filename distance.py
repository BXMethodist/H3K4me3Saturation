"""
This module is different from selector which choose the best match based on overlap
This module is based on the distance, try to look for the distance for each member
"""
from bisect import bisect
import os, numpy as np, pandas as pd, pickle
from converter import gene
from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit


def gene_to_region_tss(feature_map, reference_map, reference_index_map, gene_annotation_table):
    """
    Map each gene to the nearest region, give them the best distance and separate the upstream and downstream
    :param feature_map:
    :param reference_map:
    :param gene_annotation_table: original annotation txt download from UCSC
    :param tss: distance to tss or distance to gene_body, if true, is distance to tss
    :return:
    """
    with open(reference_map, 'rb') as f:
        reference_map = pickle.load(f)
    f.close()

    with open(reference_index_map, 'rb') as f:
        reference_index_map = pickle.load(f)
    f.close()

    with open(feature_map, 'rb') as f:
        featuremap = pickle.load(f)
    f.close()

    upstream_results = []
    downstream_results = []

    gtf_df = pd.read_csv(gene_annotation_table, sep='\t')

    count = 0

    for gene_name2 in gtf_df['name2'].unique():
        cur_df = gtf_df[gtf_df['name2'] == gene_name2]
        chromosome = cur_df.iloc[0]['chrom']
        if chromosome not in reference_index_map:
            count += 1
            continue
        indexes = reference_index_map[chromosome]
        best_distance = None
        up_down = None
        for i in range(cur_df.shape[0]):
            cur_tss = cur_df.iloc[i]['txStart']
            insert_indexes = bisect(indexes, cur_tss)
            left = None
            right = None
            if insert_indexes -1 >=0:
                left = indexes[insert_indexes -1]
            if insert_indexes <= len(indexes) -1:
                right = indexes[insert_indexes]
            candidate_regions = set()
            if left is not None:
                candidate_regions = candidate_regions.union(reference_map[chromosome][left])
            if right is not None:
                candidate_regions = candidate_regions.union(reference_map[chromosome][right])
            for region in candidate_regions:
                if region.start <= cur_tss <= region.end:
                    distance = 0
                else:
                    if abs(region.start - cur_tss) >= abs(region.end - cur_tss):
                        distance = region.end - cur_tss
                    else:
                        distance = region.start - cur_tss

                # mid = (region.end + region.start)/2
                # distance = mid - cur_tss
                if best_distance is None or abs(distance) < best_distance:
                    best_distance = abs(distance)
                    if cur_df.iloc[i]['strand'] == '+':
                        if distance <= 0:
                            up_down = -1
                        else:
                            up_down = 1
                    if cur_df.iloc[i]['strand'] == '-':
                        if distance >= 0:
                            up_down = -1
                        else:
                            up_down = 1
        if up_down == -1:
            upstream_results.append((gene_name2, best_distance))
        elif up_down == 1:
            downstream_results.append((gene_name2, best_distance))

    total_results = upstream_results + downstream_results

    total_results = sorted(total_results, key=lambda x:x[1])
    upstream_results = sorted(upstream_results, key=lambda x: x[1])
    downstream_results = sorted(downstream_results, key=lambda x: x[1])

    print count

    total_results_df = pd.DataFrame(total_results)
    total_results_df.to_csv('total_gene_distance.csv', index=False, header=False)

    upstream_results = pd.DataFrame(upstream_results)
    upstream_results.to_csv('upstream_gene_distance.csv', index=False, header=False)

    downstream_results = pd.DataFrame(downstream_results)
    downstream_results.to_csv('downstream_gene_distance.csv', index=False, header=False)

    return total_results, upstream_results, downstream_results


def transcript_to_variant_tss(feature_map, variant_map, variant_index_map, gene_annotation_table):
    """
    Map each transcript to the nearest variant, give them the best distance and separate the upstream and downstream
    :param feature_map:
    :param reference_map:
    :return:
    """
    with open(variant_map, 'rb') as f:
        variant_map = pickle.load(f)
    f.close()

    with open(variant_index_map, 'rb') as f:
        variant_index_map = pickle.load(f)
    f.close()

    with open(feature_map, 'rb') as f:
        featuremap = pickle.load(f)
    f.close()

    upstream_results = []
    downstream_results = []

    gtf_df = pd.read_csv(gene_annotation_table, sep='\t')

    count = 0

    for gene_name in gtf_df['name'].unique():
        cur_df = gtf_df[gtf_df['name'] == gene_name]
        chromosome = cur_df.iloc[0]['chrom']
        if chromosome not in variant_index_map:
            count += 1
            continue
        indexes = variant_index_map[chromosome]
        best_distance = None
        up_down = None
        for i in range(cur_df.shape[0]):
            cur_tss = cur_df.iloc[i]['txStart']
            insert_indexes = bisect(indexes, cur_tss)
            left = None
            right = None
            if insert_indexes -1 >=0:
                left = indexes[insert_indexes -1]
            if insert_indexes<= len(indexes) -1:
                right = indexes[insert_indexes]
            candidate_variants = set()
            if left is not None:
                candidate_variants = candidate_variants.union(variant_map[chromosome][left])
            if right is not None:
                candidate_variants = candidate_variants.union(variant_map[chromosome][right])
            for variant in candidate_variants:
                if variant.left_boundary <= cur_tss <= variant.right_boundary:
                    distance = 0
                else:
                    if abs(variant.left_boundary - cur_tss) >= abs(variant.right_boundary - cur_tss):
                        distance = variant.right_boundary - cur_tss
                    else:
                        distance = variant.left_boundary - cur_tss

                # mid = variant.center
                # distance = mid - cur_tss
                if best_distance is None or abs(distance) < best_distance:
                    best_distance = abs(distance)
                    if cur_df.iloc[i]['strand'] == '+':
                        if distance <= 0:
                            up_down = -1
                        else:
                            up_down = 1
                    if cur_df.iloc[i]['strand'] == '-':
                        if distance >= 0:
                            up_down = -1
                        else:
                            up_down = 1
        if up_down == -1:
            upstream_results.append((gene_name, best_distance))
        elif up_down == 1:
            downstream_results.append((gene_name, best_distance))

    total_results = upstream_results + downstream_results

    total_results = sorted(total_results, key=lambda x:x[1])
    upstream_results = sorted(upstream_results, key=lambda x: x[1])
    downstream_results = sorted(downstream_results, key=lambda x: x[1])

    print count

    total_results_df = pd.DataFrame(total_results)
    total_results_df.to_csv('total_transcript_distance.csv', index=False, header=False)

    upstream_results = pd.DataFrame(upstream_results)
    upstream_results.to_csv('upstream_transcript_distance.csv', index=False, header=False)

    downstream_results = pd.DataFrame(downstream_results)
    downstream_results.to_csv('downstream_transcript_distance.csv', index=False, header=False)

    return total_results, upstream_results, downstream_results

def region_to_gene(feature_map, feature_index_map, referencelistmap):
    """
    Map each variant to the nearest transcript, give them the best distance and separate the upstream and downstream
    :param feature_map: feature map with tss only
    :param feature_index_map: feature index map with tss only
    :param referencelistmap: the reference map contains the list of region which is the original
    reference map generated by RefRegion.
    :return:
    """
    with open(feature_index_map, 'rb') as f:
        feature_index_map = pickle.load(f)
    f.close()

    with open(feature_map, 'rb') as f:
        featuremap = pickle.load(f)
    f.close()

    with open(referencelistmap, 'rb') as f:
        referencelistmap = pickle.load(f)
    f.close()

    upstream_results = []
    downstream_results = []

    count = 0

    for region in referencelistmap:
        cur_left_boundary = region.start
        cur_right_boundary = region.end
        cur_chr = region.chromosome

        if cur_chr not in feature_index_map:
            count += 1
            continue
        else:
            indexes = feature_index_map[cur_chr]

        left_index = bisect(indexes, cur_left_boundary)
        right_index = bisect(indexes, cur_right_boundary)

        # the rational is if overlap, the distance is zero
        if left_index < right_index:
            upstream_results.append((region.id, 0))
            continue
        else:
            potential_indexes = []
            if left_index - 1>=0:
                potential_indexes.append(indexes[left_index-1])
            if right_index <= len(indexes) -1:
                potential_indexes.append(indexes[right_index])
            candidate_features = set()
            for feature_position in potential_indexes:
                candidate_features = candidate_features.union(feature_map[cur_chr][feature_position])

        best_distance = None
        up_down = None
        for feature in candidate_features:
            cur_tss = feature.tss
            if region.start <= cur_tss <= region.end:
                distance = 0
            else:
                if abs(region.start - cur_tss) >= abs(region.end - cur_tss):
                    distance = region.end - cur_tss
                else:
                    distance = region.start - cur_tss

            # mid = variant.center
            # distance = mid - cur_tss
            if best_distance is None or abs(distance) < best_distance:
                best_distance = abs(distance)
                if feature.strand == '+':
                    if distance <= 0:
                        up_down = -1
                    else:
                        up_down = 1
                if feature.strand == '-':
                    if distance >= 0:
                        up_down = -1
                    else:
                        up_down = 1
        if up_down == -1:
            upstream_results.append((region.id, best_distance))
        elif up_down == 1:
            downstream_results.append((region.id, best_distance))

    total_results = upstream_results + downstream_results

    total_results = sorted(total_results, key=lambda x:x[1])
    upstream_results = sorted(upstream_results, key=lambda x: x[1])
    downstream_results = sorted(downstream_results, key=lambda x: x[1])

    print count

    total_results_df = pd.DataFrame(total_results)
    total_results_df.to_csv('total_region_distance.csv', index=False, header=False)

    upstream_results = pd.DataFrame(upstream_results)
    upstream_results.to_csv('upstream_region_distance.csv', index=False, header=False)

    downstream_results = pd.DataFrame(downstream_results)
    downstream_results.to_csv('downstream_region_distance.csv', index=False, header=False)

    return total_results, upstream_results, downstream_results

def variant_to_transcript_tss(feature_map, feature_index_map, referencelistmap):
    """
    Map each variant to the nearest transcript, give them the best distance and separate the upstream and downstream
    :param feature_map: feature map with tss only
    :param feature_index_map: feature index map with tss only
    :param referencelistmap: the reference map contains the list of region which is the original
    reference map generated by RefRegion.
    :return:
    """
    with open(feature_index_map, 'rb') as f:
        feature_index_map = pickle.load(f)
    f.close()

    with open(feature_map, 'rb') as f:
        featuremap = pickle.load(f)
    f.close()

    with open(referencelistmap, 'rb') as f:
        referencelistmap = pickle.load(f)
    f.close()

    upstream_results = []
    downstream_results = []

    count = 0

    for region in referencelistmap:
        for variant in region.variants:
            cur_left_boundary = variant.left_boundary
            cur_right_boundary = variant.right_boundary
            cur_chr = variant.chromosome

            if cur_chr not in feature_index_map:
                count += 1
                continue
            else:
                indexes = feature_index_map[cur_chr]

            left_index = bisect(indexes, cur_left_boundary)
            right_index = bisect(indexes, cur_right_boundary)

            # the rational is if overlap, the distance is zero
            if left_index < right_index:
                upstream_results.append((variant.id, 0))
                continue
            else:
                potential_indexes = []
                if left_index - 1>=0:
                    potential_indexes.append(indexes[left_index-1])
                if right_index <= len(indexes) -1:
                    potential_indexes.append(indexes[right_index])
                candidate_features = set()
                for feature_position in potential_indexes:
                    candidate_features = candidate_features.union(feature_map[cur_chr][feature_position])

            best_distance = None
            up_down = None
            for feature in candidate_features:
                cur_tss = feature.tss
                if variant.left_boundary <= cur_tss <= variant.right_boundary:
                    distance = 0
                else:
                    if abs(variant.left_boundary - cur_tss) >= abs(variant.right_boundary - cur_tss):
                        distance = variant.right_boundary - cur_tss
                    else:
                        distance = variant.left_boundary - cur_tss

                # mid = variant.center
                # distance = mid - cur_tss
                if best_distance is None or abs(distance) < best_distance:
                    best_distance = abs(distance)
                    if feature.strand == '+':
                        if distance <= 0:
                            up_down = -1
                        else:
                            up_down = 1
                    if feature.strand == '-':
                        if distance >= 0:
                            up_down = -1
                        else:
                            up_down = 1
            if up_down == -1:
                upstream_results.append((variant.id, best_distance))
            elif up_down == 1:
                downstream_results.append((variant.id, best_distance))

    total_results = upstream_results + downstream_results

    total_results = sorted(total_results, key=lambda x:x[1])
    upstream_results = sorted(upstream_results, key=lambda x: x[1])
    downstream_results = sorted(downstream_results, key=lambda x: x[1])

    print count

    total_results_df = pd.DataFrame(total_results)
    total_results_df.to_csv('total_variant_distance.csv', index=False, header=False)

    upstream_results = pd.DataFrame(upstream_results)
    upstream_results.to_csv('upstream_variant_distance.csv', index=False, header=False)

    downstream_results = pd.DataFrame(downstream_results)
    downstream_results.to_csv('downstream_variant_distance.csv', index=False, header=False)

    return total_results, upstream_results, downstream_results


# gene_to_region('./pkl/hg19_RefSeq_refGenetss_only.pkl', './pkl/75_distance_3kb.pkl', './pkl/75_distance_3kb_index.pkl',
#                './pkl/hg19_RefSeq_refGene.txt')




