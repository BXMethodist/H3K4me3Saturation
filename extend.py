"""
This module is used to extend a reference map based on the other reference map
"""

import pickle, pandas as pd
from bisect import bisect
from collections import defaultdict
from clusterUtils import genome_size_chrom

def isOverlap(peak, ref_distance_map, ref_distance_indexmap):
    """
    check whether this peak need to be keep or not
    :param peak: the peak from danpos file,
    :param ref_distance_map: the refmap object
    :param ref_distance_indexmap: the refmap object
    :return:
    """
    chromosome = peak[0]
    start = int(peak[1])
    end = int(peak[2])

    if chromosome not in ref_distance_indexmap:
        return False

    indexes = ref_distance_indexmap[chromosome]

    left_index = bisect(indexes, start)
    right_index = bisect(indexes, end)

    # the rational is if overlap, the distance is zero
    candidate_regions = set()

    potential_indexes = []

    left_index = left_index - 10 if left_index - 10 >= 0 else 0
    for index in indexes[left_index - 1: right_index+10]:
        potential_indexes.append(index)

    for feature_position in potential_indexes:
        candidate_regions = candidate_regions.union(ref_distance_map[chromosome][feature_position])

    for region in candidate_regions:
        if start <= region.start <= end:
            return True
        if start <= region.end <= end:
            return True
        if region.start <= start and end <= region.end:
            return True
    return False

def filter_peaks(peak_table, ref_distance_map, ref_distance_indexmap, output='./extend_100_10/'):
    """
    get a table of peaks and return a table with filtered peaks
    :param peak_table: danpos peak table
    :param ref_distance_map: refmap object
    :param ref_distance_indexmap: refmap object
    :return:
    """
    peaks = open(peak_table, 'r')
    info = peaks.readlines()
    peaks.close()

    results = []

    for line in info:
        line = line.strip()
        line = line.split('\t')
        if line[1] == "start":
            results.append('\t'.join(line))
            continue
        if isOverlap(line, ref_distance_map, ref_distance_indexmap):
            results.append('\t'.join(line))

    file_name = peak_table[peak_table.rfind('/')+1:]

    f = open(output+file_name, "w")
    for line in results:
        f.write(line+'\n')
    f.close()

def extend_peaks_danpos(seed_peak_table_path, seed, target, cutoff=2200):
    """
    :param seed_peak_table: danpos peak table serve as seed
    :param extend_peak_table: danpos peak table serve as extension
    :return:
    """
    seed_peak_table = open(seed_peak_table_path, 'r')
    info1 = seed_peak_table.readlines()
    seed_peak_table.close()

    seed_peaks = {}

    for line1 in info1:
        line1 = line1.strip()
        line1 = line1.split('\t')
        if line1[1] == "start":
            continue

        chromosome = line1[0]
        if chromosome not in seed_peaks:
            seed_peaks[chromosome] = []
        seed_peaks[chromosome].append(line1)

    extend_peak_tables = {}

    for j in range(target, seed, 10):
        extend_peak_table_path = seed_peak_table_path.replace(str(seed)+'cutoff', str(j)+'cutoff')
        extend_peak_table = open(extend_peak_table_path, 'r')
        info2 = extend_peak_table.readlines()
        extend_peak_table.close()

        extend_peaks = {}

        for line2 in info2:
            line2 = line2.strip()
            line2 = line2.split('\t')
            if line2[1] == "start":
                continue

            chromosome = line2[0]
            if chromosome not in extend_peaks:
                extend_peaks[chromosome] = []
            extend_peaks[chromosome].append(line2)
        for chromosome in seed_peaks.keys():
            if chromosome in extend_peaks:
                extend_peaks[chromosome] = sorted(extend_peaks[chromosome], key=lambda y: int(y[1]))
        extend_peak_tables[j] = extend_peaks

    for chromosome in seed_peaks.keys():
        seed_peaks[chromosome] = sorted(seed_peaks[chromosome], key=lambda x: int(x[1]))

    result_extensions = set()
    distance_extensions = []
    percentage = []
    overlapped = 0
    not_overlapped = 0
    extend_peaks = extend_peak_tables[target]
    for chromosome in seed_peaks.keys():
        if chromosome in extend_peaks:
            bisect_keys = [int(z[1]) for z in extend_peaks[chromosome]]
            for peak in seed_peaks[chromosome]:
                start, end = peak[1], peak[2]
                index = bisect(bisect_keys, int(start))
                if index - 1 >=0:
                    left = index - 1
                else:
                    left = 0
                if index < len(bisect_keys):
                    right = index
                else:
                    right = len(bisect_keys) - 1
                added = False
                candidates = [extend_peaks[chromosome][left], extend_peaks[chromosome][right]]
                for candidate in candidates:
                    cur_start, cur_end = candidate[1], candidate[2]
                    if int(cur_start) <= int(start) <= int(end) <= int(cur_end):
                        if int(start) - int(cur_start) + int(cur_end) - int(end) > cutoff:
                            new_peak = find_right_peak(extend_peak_tables, chromosome, start, end, cutoff)
                            if new_peak is not None:
                                result_extensions.add(tuple(new_peak))
                            else:
                                result_extensions.add(tuple(peak))
                        else:
                            result_extensions.add(tuple(candidate))
                        # distance_extensions.append(int(start)-int(cur_start)+int(cur_end)-int(end))
                        # percentage.append((int(start)-int(cur_start)+int(cur_end)-int(end))*1.0/(int(end)-int(start)))
                        added = True
                        break
                if not added:
                    print candidates
                    print peak
                    not_overlapped += 1
                else:
                    overlapped += 1
        else:
            for peak in seed_peaks[chromosome]:
                result_extensions.add(tuple(peak))


    file_name = seed_peak_table_path[seed_peak_table_path.rfind('/')+1:]
    extend_file = open('./extend_100_10/'+file_name, 'w')
    for line in result_extensions:
        # line = [str(l) for l in line]
        line = '\t'.join(line) + '\n'
        extend_file.write(line)
    extend_file.close()

    print not_overlapped, overlapped
    return percentage

def compare_peak_number(extended_path, original_path):
    extend_tables = [x for x in os.listdir(extended_path)]
    results = []
    for i in range(len(extend_tables)):
        extend_table = extended_path + '/' + extend_tables[i]
        original_table = original_path + '/' + extend_tables[i]

        e = open(extend_table, 'r')
        info_e = e.readlines()
        if len(info_e) == 0:
            print extend_table
            continue
        e.close()
        first_e = info_e[0].split('\t')[1]
        if first_e == 'start':
            len_e = len(info_e) - 1
        else:
            len_e = len(info_e)

        o = open(original_table, 'r')
        info_o = o.readlines()
        o.close()
        if len(info_o) == 0:
            print original_table
            continue
        first_o = info_o[0].split('\t')[1]
        if first_o == 'start':
            len_o = len(info_o) - 1
        else:
            len_o = len(info_o)
        results.append((len_e, len_o))

    df = pd.DataFrame(results)
    df.to_csv('extend_peak_number_comp_100_10_vs_100.csv', index=False, header=None)
    return df

def extend_width_distribution(sample):
    f = open(sample, 'r')
    info = f.readlines()
    f.close()

    results = []

    for line1 in info:
        line1 = line1.strip()
        line1 = line1.split('\t')
        if line1[1] == "start":
            continue
        start = int(line1[1])
        end = int(line1[2])
        results.append([end-start])

    return results

def get_all_peaks(path):
    tables = os.listdir(path)
    peaks = defaultdict(list)
    for table in tables:
        f = open(path+'/'+table, 'r')
        info = f.readlines()
        f.close()
        for line in info:
            line = line.strip().split('\t')
            if line[1] == 'start':
                continue
            else:
                peaks[line[0]].append(line[1:]+[table])
    for key in peaks.keys():
        peaks[key] = sorted(peaks[key], key=lambda x:int(x[0]))
    return peaks

def get_peaks(chromosome, start, end, peaks):
    """
    get all peaks in this region
    :param chromosome:
    :param start:
    :param end:
    :return:
    """
    candidates = peaks[chromosome]
    indexes = [int(x[0]) for x in candidates]
    index_start = bisect(indexes, start)
    index_end = bisect(indexes, end)
    return candidates[index_start:index_end]

def region_peaks_number(refmap, peaks):
    input_file = open(refmap, "r")
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

    regular_chrs = genome_size_chrom()

    results = []

    for key in regular_chrs.keys():
        for peak in region_map[key]:
            start, end = peak
            results.append(len(get_peaks(key, start, end, peaks)))

    f = open(refmap[:-4]+'_regions_per_samples.txt', 'w')
    for line in results:
        f.write(str(line) + '\n')
    f.close()

def find_right_peak(extend_peak_tables, chromosome, start, end, cutoff):
    """
    :param seed_peak_table_path:
    :param chromosome:
    :param start:
    :param end:
    :param cutoff:
    :return:
    """
    for i in range(20, 100, 10):
        seed_peaks = extend_peak_tables[i]
        bisect_keys = [int(z[1]) for z in seed_peaks[chromosome]]

        index = bisect(bisect_keys, int(start))
        if index - 1 >= 0:
            left = index - 1
        else:
            left = 0
        if index < len(bisect_keys):
            right = index
        else:
            right = len(bisect_keys) - 1
        candidates = [seed_peaks[chromosome][left], seed_peaks[chromosome][right]]
        for candidate in candidates:
            cur_start, cur_end = candidate[1], candidate[2]
            if int(cur_start) <= int(start) <= int(end) <= int(cur_end):
                if int(start) - int(cur_start) + int(cur_end) - int(end) > cutoff:
                    continue
                else:
                    return candidate
    return None




import os

#
seed_peak_tables = ['/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/'+x for x in os.listdir('/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/')]
#
#
#
distance_extensions = []

for i in range(len(seed_peak_tables)):
    seed_peak_table = seed_peak_tables[i]
    extend_peak_table = seed_peak_table.replace('100cutoff', '10cutoff')
    print seed_peak_table
    print extend_peak_table
    distance_extensions += extend_peaks_danpos(seed_peak_table, 100, 10)




# extend_distance_file = open('./extend_distance_percentage.txt', 'w')
# for line in distance_extensions:
#     # line = [str(l) for l in line]
#     extend_distance_file.write(str(line)+'\n')
# extend_distance_file.close()

# compare_peak_number("/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10", "/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/")

# seed_peak_tables = ['/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/'+x for x in os.listdir('/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10')]
#
# extend_results = []
# for table in seed_peak_tables:
#     extend_results += extend_width_distribution(table)
#
# df = pd.DataFrame(extend_results)
# df.to_csv('extend_raw_peak.csv', index=False, header=None)
#
# seed_peak_tables = ['/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/'+x for x in os.listdir('/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/')]
#
# original_results = []
# for table in seed_peak_tables:
#     original_results += extend_width_distribution(table)
#
# df = pd.DataFrame(original_results)
# df.to_csv('original_raw_peak.csv', index=False, header=None)

# cutoffs = [x for x in range(10, 100,10)]
#
# for cutoff in cutoffs:
#     peaks = get_all_peaks('/home/tmhbxx3/archive/KFH3K4me3/'+str(cutoff)+'cutoff/pooled')
#     refmap = str(cutoff) +"_refmap.csv"
#     region_peaks_number(refmap, peaks)







