import os
import collections

path = "/Users/boxia/Desktop/pictures/"
output = "/Users/boxia/Desktop/single_clusters/"

list_of_pictures = [x[:-13] for x in os.listdir(path)]

single_pictures = [item for item, count in collections.Counter(list_of_pictures).items() if count == 1]

for picture in single_pictures:
    cmd = "mv "+path+picture+"* "+output
    os.system(cmd)
