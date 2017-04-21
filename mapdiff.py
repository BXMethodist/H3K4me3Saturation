"""
This module is used to test the difference between the reference map
"""

from RefRegion import ReferenceRegion, ReferenceVariant, ReferenceUnit
from simple_region import SimpleRegion
import os, pickle, bisect
from collections import defaultdict


def difference_maps(large_distance_map, small_distance_map):
    """

    :param large_distance_map:
    :param small_distance_map:
    :return: the lists of group which the region are merged
    """
    overlapped = defaultdict(set)
    for key in small_distance_map.keys():
        if key not in large_distance_map:
            continue
        indexes = sorted(large_distance_map[key].keys())
        for value in small_distance_map[key].keys():
            insert_index = bisect.bisect(indexes, value)
            if insert_index - 1 >= 0:
                insert_index = insert_index - 1
            else:
                insert_index = 0
            cur_pos = indexes[insert_index]
            cur_large_region = list(large_distance_map[key][cur_pos])[0]
            if len(large_distance_map[key][cur_pos]) > 1:
                print large_distance_map[key][cur_pos]
            overlapped[(cur_large_region.chromosome, cur_large_region.start, cur_large_region.end)].add(list(small_distance_map[key][value])[0])

    for key in overlapped.keys():
        if len(overlapped[key]) == 1:
            del overlapped[key]
    # print len(overlapped.keys())

    merge_distance = []
    for value in overlapped.values():
        cur_value = sorted(list(value), key=lambda x : x.start)
        for i in range(len(cur_value) - 1):
            cur = cur_value[i]
            next_cur = cur_value[i+1]
            merge_distance.append(next_cur.start - cur.end)
    return merge_distance







# with open('./ref100/100_distance.pkl', 'rb') as f:
#     small_map = pickle.load(f)
# f.close()
#
# with open('./extend_100_10/extend_100_10_simple_region_distance.pkl', 'rb') as f:
#     large_map = pickle.load(f)
# f.close()
#
#
# merge_dis = difference_maps(large_map, small_map)

# with open('100_10_extend_merge_distance' + '.pkl', 'wb') as f:
#     pickle.dump(merge_dis, f, pickle.HIGHEST_PROTOCOL)
# f.close()

with open('100_10_extend_merge_distance.pkl', 'rb') as f:
    merge_dis = pickle.load(f)
f.close()

merge_dis = sorted(merge_dis)

f = open('100_10_extend_merge_distance.txt', 'w')
for line in merge_dis:
    f.write(str(line)+'\n')
f.close()

# import matplotlib.pyplot as plt
#
# plt.hist(merge_dis, bins=100)
# plt.show()


