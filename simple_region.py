"""
This module is a simply region object, similar as the one generated from refregion, but it only has start, end and chromosome.
"""

import os, pickle
from collections import defaultdict

def get_region_map(map):
    map_obj = open(map, 'r')
    info = map_obj.readlines()
    map_obj.close()

    regions = []
    for line in info:
        if line.startswith('>'):
            chromosome = line.strip()[1:]
        else:
            start, end = line.strip().split(',')
            cur_region = SimpleRegion(chromosome, int(start)*10, int(end)*10)
            regions.append(cur_region)

    with open('extend_100_10_simple_region' + '.pkl', 'wb') as f:
        pickle.dump(regions, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    print len(regions)
    regions_distance = {}
    regions_index = defaultdict(set)
    for region in regions:
        if region.chromosome not in regions_distance:
            regions_distance[region.chromosome] = defaultdict(set)
        regions_distance[region.chromosome][region.start].add(region)
        regions_index[region.chromosome].add(region.start)
    # print regions_distance
    for key in regions_index.keys():
        regions_index[key] = sorted(list(regions_index[key]))
    with open('extend_100_10_2200_simple_region_distance' + '.pkl', 'wb') as f:
        pickle.dump(regions_distance, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    with open('extend_100_10_2200_simple_region_distance_index' + '.pkl', 'wb') as f:
        pickle.dump(regions_index, f, pickle.HIGHEST_PROTOCOL)
    f.close()
    return

class SimpleRegion():
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.id = chromosome[3:] + '.' + str(start) +'.' +str(end)

# get_region_map('/Users/boxia/PycharmProjects/H3K4me3Saturation/extend_100_10_2200/100_10_2200_refmap.csv')








