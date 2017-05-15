"""
This module is different from selector which choose the best match based on overlap
This module is based on the distance, try to look for the distance for each member
"""
from bisect import bisect
import os, numpy as np, pandas as pd, pickle
from converter import gene
from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit
from simple_region import SimpleRegion
from collections import defaultdict
from itertools import combinations


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

    results = []

    gtf_df = pd.read_csv(gene_annotation_table, sep='\t')

    # count = 0
    # print reference_index_map.keys()

    for gene_name2 in gtf_df['name2'].unique():
        # print gene_name2
        cur_df = gtf_df[gtf_df['name2'] == gene_name2]
        chromosome = cur_df.iloc[0]['chrom']
        if chromosome not in reference_index_map:
            # count += 1
            continue
        indexes = reference_index_map[chromosome]
        best_up_distance = None

        best_down_distance = None

        best_distance = None

        for i in range(cur_df.shape[0]):
            if cur_df.iloc[i]['strand'] == '+':
                cur_tss = int(cur_df.iloc[i]['txStart'])
            elif cur_df.iloc[i]['strand'] == '-':
                cur_tss = int(cur_df.iloc[i]['txEnd'])
            insert_indexes = bisect(indexes, cur_tss)

            left = 1
            candidate_regions = set()
            while insert_indexes - left >= 0 and left < 10:
                candidate_regions = candidate_regions.union(reference_map[chromosome][indexes[insert_indexes - left]])
                left +=1

            right = 0
            while insert_indexes + right < len(indexes) and right < 10:
                candidate_regions = candidate_regions.union(reference_map[chromosome][indexes[insert_indexes+right]])
                right += 1

            # if gene_name2 == 'NOTCH1':
            #     print [(x.start, x.end) for x in candidate_regions]

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
                if cur_df.iloc[i]['strand'] == '+':
                    if distance < 0:
                        if best_up_distance is None or abs(distance) < best_up_distance:
                            best_up_distance = abs(distance)
                    elif distance > 0:
                        if best_down_distance is None or abs(distance) < best_down_distance:
                            best_down_distance = abs(distance)
                    else:
                        best_up_distance = distance
                        best_down_distance = distance
                if cur_df.iloc[i]['strand'] == '-':
                    if distance > 0:
                        if best_up_distance is None or abs(distance) < best_up_distance:
                            best_up_distance = abs(distance)
                    elif distance < 0:
                        if best_down_distance is None or abs(distance) < best_down_distance:
                            best_down_distance = abs(distance)
                    else:
                        best_up_distance = distance
                        best_down_distance = distance
                if best_distance is None or abs(distance) < best_distance:
                    best_distance = abs(distance)
        # print best_down_distance, best_up_distance, len(upstream_results), len(downstream_results)
        results.append((gene_name2, best_up_distance, best_down_distance, best_distance))

        if best_up_distance is None or best_down_distance is None:
            print 'something wrong!'

    # print count

    results_df = pd.DataFrame(results)
    results_df.to_csv('100_10_gene_distance.csv', index=False, header=False)

    return results_df

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

    # count = 0

    for gene_name in gtf_df['name'].unique():
        cur_df = gtf_df[gtf_df['name'] == gene_name]
        chromosome = cur_df.iloc[0]['chrom']
        if chromosome not in variant_index_map:
            # count += 1
            continue
        indexes = variant_index_map[chromosome]
        best_up_distance = None
        best_down_distance = None

        for i in range(cur_df.shape[0]):
            if cur_df.iloc[i]['strand'] == '+':
                cur_tss = int(cur_df.iloc[i]['txStart'])
            elif cur_df.iloc[i]['strand'] == '-':
                cur_tss = int(cur_df.iloc[i]['txEnd'])
            insert_indexes = bisect(indexes, cur_tss)

            left = 1
            candidate_variants = set()
            while insert_indexes - left >= 0 and left < 10:
                candidate_variants = candidate_variants.union(variant_map[chromosome][indexes[insert_indexes - left]])
                left += 1

            right = 0
            while insert_indexes + right < len(indexes) and right < 10:
                candidate_variants = candidate_variants.union(variant_map[chromosome][indexes[insert_indexes + right]])
                right += 1

            for variant in candidate_variants:
                if variant.left_boundary <= cur_tss <= variant.right_boundary:
                    distance = 0
                else:
                    if abs(variant.left_boundary - cur_tss) >= abs(variant.right_boundary - cur_tss):
                        distance = variant.right_boundary - cur_tss
                    else:
                        distance = variant.left_boundary - cur_tss

                if cur_df.iloc[i]['strand'] == '+':
                    if distance <= 0:
                        if best_up_distance is None or abs(distance) < best_up_distance:
                            best_up_distance = abs(distance)
                    else:
                        if best_down_distance is None or abs(distance) < best_down_distance:
                            best_down_distance = abs(distance)
                if cur_df.iloc[i]['strand'] == '-':
                    if distance >= 0:
                        if best_up_distance is None or abs(distance) < best_up_distance:
                            best_up_distance = abs(distance)
                    else:
                        if best_down_distance is None or abs(distance) < best_down_distance:
                            best_down_distance = abs(distance)

        upstream_results.append((gene_name, best_up_distance))
        downstream_results.append((gene_name, best_down_distance))

    upstream_results = sorted(upstream_results, key=lambda x: x[1])
    downstream_results = sorted(downstream_results, key=lambda x: x[1])

    # print count
    upstream_results = pd.DataFrame(upstream_results)
    upstream_results.to_csv('upstream_transcript_distance.csv', index=False, header=False)

    downstream_results = pd.DataFrame(downstream_results)
    downstream_results.to_csv('downstream_transcript_distance.csv', index=False, header=False)

    return upstream_results, downstream_results

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
        feature_map = pickle.load(f)
    f.close()

    with open(referencelistmap, 'rb') as f:
        referencelistmap = pickle.load(f)
    f.close()

    results = []

    count = 0

    # print len(referencelistmap)
    # print feature_index_map.keys()

    for region in referencelistmap:
        cur_left_boundary = region.start
        cur_right_boundary = region.end
        cur_chr = region.chromosome

        # print cur_chr

        if cur_chr not in feature_index_map:
            count += 1
            continue
        else:
            indexes = feature_index_map[cur_chr]

        left_index = bisect(indexes, cur_left_boundary)
        right_index = bisect(indexes, cur_right_boundary)

        # the rational is if overlap, the distance is zero
        candidate_features = set()

        potential_indexes = []

        left_index = left_index -20 if left_index -20 >=0 else 0
        for index in indexes[left_index: right_index+20]:
            potential_indexes.append(index)

        for feature_position in potential_indexes:
            candidate_features = candidate_features.union(feature_map[cur_chr][feature_position])

        best_up_distance = None
        best_down_distance = None
        best_distance = None
        distance = None

        for feature in candidate_features:
            cur_tss = feature.tss
            if region.start <= cur_tss <= region.end:
                distance = 0
            else:
                if abs(region.start - cur_tss) >= abs(region.end - cur_tss):
                    distance = region.end - cur_tss
                else:
                    distance = region.start - cur_tss

            if best_distance is None or abs(distance) < best_distance:
                best_distance = abs(distance)

            if feature.strand == 1:
                if distance < 0:
                    if best_up_distance is None or abs(distance) < best_up_distance:
                        best_up_distance = abs(distance)
                elif distance > 0:
                    if best_down_distance is None or abs(distance) < best_down_distance:
                        best_down_distance = abs(distance)
                else:
                    # print "why?", distance ==0
                    best_down_distance = distance
                    best_up_distance = distance
            elif feature.strand == -1:
                if distance > 0:
                    if best_up_distance is None or abs(distance) < best_up_distance:
                        best_up_distance = abs(distance)
                elif distance < 0:
                    if best_down_distance is None or abs(distance) < best_down_distance:
                        best_down_distance = abs(distance)
                else:
                    best_up_distance = distance
                    best_down_distance = distance
            else:
                print feature.strand, 'feature'
            # print best_up_distance, best_down_distance
        # print region.id, best_distance
        results.append((region.id, best_up_distance, best_down_distance, best_distance))
        if best_up_distance is None or best_down_distance is None or best_up_distance < 0 and best_down_distance<0:
            print "something wrong!"

    # print len(upstream_results), len(downstream_results)

    results_df = pd.DataFrame(results)
    results_df.to_csv('100_10_region_distance.csv', index=False, header=False)

    return results_df

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
        feature_map = pickle.load(f)
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
            candidate_features = set()

            potential_indexes = []

            left_index = left_index - 20 if left_index - 20 >= 0 else 0
            for index in indexes[left_index: right_index+20]:
                potential_indexes.append(index)

            for feature_position in potential_indexes:
                candidate_features = candidate_features.union(feature_map[cur_chr][feature_position])

            best_up_distance = None
            best_down_distance = None

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
                if feature.strand == 1:
                    if distance <= 0:
                        if best_up_distance is None or abs(distance) < best_up_distance:
                            best_up_distance = abs(distance)
                    else:
                        if best_down_distance is None or abs(distance) < best_down_distance:
                            best_down_distance = abs(distance)
                if feature.strand == -1:
                    if distance >= 0:
                        if best_up_distance is None or abs(distance) < best_up_distance:
                            best_up_distance = abs(distance)
                    else:
                        if best_down_distance is None or abs(distance) < best_down_distance:
                            best_down_distance = abs(distance)

            upstream_results.append((variant.id, best_up_distance))
            downstream_results.append((variant.id, best_down_distance))

    upstream_results = sorted(upstream_results, key=lambda x: x[1])
    downstream_results = sorted(downstream_results, key=lambda x: x[1])

    # print count

    upstream_results = pd.DataFrame(upstream_results)
    upstream_results.to_csv('upstream_variant_distance.csv', index=False, header=False)

    downstream_results = pd.DataFrame(downstream_results)
    downstream_results.to_csv('downstream_variant_distance.csv', index=False, header=False)

    return upstream_results, downstream_results

def sample_tss(feature_map, feature_index_map, peak_table):
    """
    this function will generate the sample peaks distance
    :param feature_map: feature map with tss only
    :param feature_index_map: feature index map with tss only
    :param peak_table: peak table for danpos
    :return:
    """
    with open(feature_index_map, 'rb') as f:
        feature_index_map = pickle.load(f)
    f.close()

    with open(feature_map, 'rb') as f:
        feature_map = pickle.load(f)
    f.close()

    file = open(peak_table, "rb")

    upstream_results = []
    downstream_results = []
    for line in file.readlines():
        info = line.split("\t")
        if info[1] == "start":
            continue
        cur_left_boundary = int(info[1]) / 10
        cur_right_boundary = int(info[2]) / 10
        cur_chr = info[0]

        peak_id = cur_chr+'.'+str(cur_left_boundary)+'.'+str(cur_right_boundary)

        if cur_chr not in feature_index_map:
            continue
        else:
            indexes = feature_index_map[cur_chr]

        left_index = bisect(indexes, cur_left_boundary)
        right_index = bisect(indexes, cur_right_boundary)

        # the rational is if overlap, the distance is zero
        if left_index < right_index:
            upstream_results.append((peak_id, 0))
            continue
        else:
            potential_indexes = []
            if left_index - 1 >= 0:
                potential_indexes.append(indexes[left_index - 1])
            if right_index <= len(indexes) - 1:
                potential_indexes.append(indexes[right_index])
            candidate_features = set()
            for feature_position in potential_indexes:
                candidate_features = candidate_features.union(feature_map[cur_chr][feature_position])

        best_distance = None
        up_down = None
        for feature in candidate_features:
            cur_tss = feature.tss
            if cur_left_boundary <= cur_tss <= cur_right_boundary:
                distance = 0
            else:
                if abs(cur_left_boundary - cur_tss) >= abs(cur_right_boundary - cur_tss):
                    distance = cur_right_boundary - cur_tss
                else:
                    distance = cur_left_boundary - cur_tss

            # mid = variant.center
            # distance = mid - cur_tss
            if best_distance is None or abs(distance) < best_distance:
                best_distance = abs(distance)
                if feature.strand == 1:
                    if distance <= 0:
                        up_down = -1
                    else:
                        up_down = 1
                if feature.strand == -1:
                    if distance >= 0:
                        up_down = -1
                    else:
                        up_down = 1
        if up_down == -1:
            upstream_results.append((peak_id, best_distance))
        elif up_down == 1:
            downstream_results.append((peak_id, best_distance))
    file.close()
    total_results = upstream_results + downstream_results

    total_results = sorted(total_results, key=lambda x:x[1])
    upstream_results = sorted(upstream_results, key=lambda x:x[1])
    downstream_results = sorted(downstream_results, key=lambda x:x[1])

    return total_results, upstream_results, downstream_results

def samples_distance(feature_map, feature_index_map, peak_tables_path):
    """
    this function will generate the samples peaks distance
    :param feature_map: feature map with tss only
    :param feature_index_map: feature index map with tss only
    :param peak_table: path that contains danpos tables
    :return:
    """
    tables = os.listdir(peak_tables_path)

    total_results = []
    upstream_results = []
    downstream_results = []

    max_total_length = 0
    max_up_length = 0
    max_down_length = 0

    for table in tables:
        table_path = peak_tables_path+'/'+table
        cur_total, cur_up, cur_down = sample_tss(feature_map, feature_index_map, table_path)

        total_results.append(cur_total)
        upstream_results.append(cur_up)
        downstream_results.append(cur_down)

        if max_total_length < len(cur_total):
            max_total_length = len(cur_total)
        if max_up_length < len(cur_up):
            max_up_length = len(cur_up)
        if max_down_length < len(cur_down):
            max_down_length = len(cur_down)

    new_total_results = []
    for total in total_results:
        diff = max_total_length - len(total)
        total = [x[1] for x in total]
        for i in range(diff):
            total.append(None)
        new_total_results.append(total)

    new_up_results = []
    for up in upstream_results:
        diff = max_up_length - len(up)
        up = [x[1] for x in up]
        for i in range(diff):
            up.append(None)
        new_up_results.append(up)

    new_down_results = []
    for down in downstream_results:
        diff = max_down_length - len(down)
        down = [x[1] for x in down]
        for i in range(diff):
            down.append(None)
        new_down_results.append(down)

    total_df = pd.DataFrame(new_total_results)
    up_df = pd.DataFrame(new_up_results)
    down_df = pd.DataFrame(new_down_results)

    total_mean = total_df.mean(axis=0)
    if len(total_mean) != total_df.shape[1]:
        print 'something wrong'
    up_mean = list(up_df.mean(axis=0))
    down_mean = list(down_df.mean(axis=0))

    result_df = pd.DataFrame(index=np.arange(len(total_mean)), columns=['total', 'up', 'down'])
    result_df['total'] = total_mean
    for i in range(len(total_mean)-len(up_mean)):
        up_mean.append(None)
    result_df['up'] = up_mean

    for i in range(len(total_mean)-len(down_mean)):
        down_mean.append(None)
    result_df['down'] = down_mean

    result_df.to_csv("samples_tss_distance.csv")
    return total_df, up_df, down_df

def distribution_region_distance(refmap):
    """
    get the distance distribution between the reference map
    :return:
    """
    input_file = open(refmap, "r")

    distance_distributions = []
    chr_name = None

    region_map = defaultdict(list)

    for line in input_file.readlines():
        if line.startswith(">"):
            chr_name = line.rstrip()[1:]
        else:
            line = line.rstrip().split(',')
            # print line
            start = int(line[0]) * 10
            end = int(line[1]) * 10
            region_map[chr_name].append((start, end))

    for key in region_map.keys():
        region_map[key] = sorted(region_map[key], key=lambda x:x[0])

    for value in region_map.values():
        indexes = [x[0] for x in value]
        for i in range(len(value)-1):
            cur_start = value[i][0]
            cur_end_index = bisect(indexes, cur_start + 100000)
            # cur_peaks = value[i: cur_end_index]
            cur_peaks = value[i: i+2]
            # print cur_peaks
            distance_distributions += peaks_distances(cur_peaks)

    distance_file = open('./region_gap_distance.txt', 'w')
    for line in distance_distributions:
        # line = [str(l) for l in line]
        distance_file.write(str(line)+'\n')
    distance_file.close()
    return distance_distributions

def distribution_peaks_distance(danpos_table):
    """
    get the distance distribution between peaks from all the samples
    :return:
    """
    peak_table = open(danpos_table, 'r')
    info1 = peak_table.readlines()
    peak_table.close()

    peaks = {}

    for line1 in info1:
        line1 = line1.strip()
        line1 = line1.split('\t')
        if line1[1] == "start":
            continue

        chromosome = line1[0]
        if chromosome not in peaks:
            peaks[chromosome] = []
        peaks[chromosome].append((int(line1[1]), int(line1[2])))

    for chromosome in peaks.keys():
        peaks[chromosome] = sorted(peaks[chromosome], key=lambda x: x[0])

    distance_distributions = []
    for value in peaks.values():
        indexes = [x[0] for x in value]
        for i in range(len(value)-1):
            cur_start = value[i][0]
            cur_end_index = bisect(indexes, cur_start + 100000)
            # cur_peaks = value[i: cur_end_index]
            cur_peaks = value[i: i+2]
            distance_distributions += peaks_distances(cur_peaks)

    return distance_distributions

def peaks_distances(peaks):
    """
    :return:
    """
    distances = []
    target_peak = peaks[0]
    for i in range(1, len(peaks)):
        cur_peak = peaks[i]
        distances.append(cur_peak[0] - target_peak[1])
    return distances

# distances = distribution_region_distance('./100_refmap.csv')



import os
#
# #
peak_tables = ['/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/' + x for x in
                    os.listdir('/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/')]
#
distances = []
#
for i in range(len(peak_tables)):
    peak_table = peak_tables[i]
    distances += distribution_peaks_distance(peak_table)
    print peak_table
#
distance_file = open('./100_refmap_peaks_distance.txt', 'w')
for line in distances:
    # line = [str(l) for l in line]
    distance_file.write(str(line)+'\n')
distance_file.close()





#
# gene_to_region_tss('./pkl/hg19_RefSeq_refGenetss_only.pkl', './extend_100_10/extend_100_10_simple_region_distance.pkl', './extend_100_10/extend_100_10_simple_region_distance_index.pkl',
#                './pkl/hg19_RefSeq_refGene.txt')


# transcript_to_variant_tss('./pkl/hg19_RefSeq_refGenetss_only.pkl', './ref100/100_variant_distance.pkl',
#                           './ref100/100_variant_distance_index.pkl',
#                           './pkl/hg19_RefSeq_refGene.txt')


# region_to_gene('./pkl/hg19_RefSeq_refGenetss_only.pkl', './pkl/hg19_RefSeq_refGenetss_only_index.pkl',
#                './extend_100_10/extend_100_10_simple_region.pkl')

# variant_to_transcript_tss('./pkl/hg19_RefSeq_refGenetss_only.pkl', './pkl/hg19_RefSeq_refGenetss_only_index.pkl',
#                './ref100/100_refregion.pkl')

# samples_distance('./pkl/hg19_RefSeq_refGenetss_only.pkl', './pkl/hg19_RefSeq_refGenetss_only_index.pkl', '/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled')