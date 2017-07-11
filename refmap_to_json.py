import json, pickle
from RefRegion import *

with open('./pkl/100_10_22000_refregions.pkl', 'rb') as f:
    regions = pickle.load(f)
f.close()

data = {}

for region in regions:
    if region.chromosome in data:
        if region.plot:
            data[region.chromosome][int(region.start)] = (int(region.start), int(region.end), 1)
            data[region.chromosome][int(region.end)] = (int(region.start), int(region.end), 1)
        else:
            data[region.chromosome][int(region.start)] = (int(region.start), int(region.end), 0)
            data[region.chromosome][int(region.end)] = (int(region.start), int(region.end), 0)
    else:
        data[region.chromosome] = {}
        if region.plot:
            data[region.chromosome][int(region.start)] = (int(region.start), int(region.end), 1)
            data[region.chromosome][int(region.end)] = (int(region.start), int(region.end), 1)
        else:
            data[region.chromosome][int(region.start)] = (int(region.start), int(region.end), 0)
            data[region.chromosome][int(region.end)] = (int(region.start), int(region.end), 0)

with open('100_10_2200_web_map.txt', 'w') as outfile:
    json.dump(data, outfile)