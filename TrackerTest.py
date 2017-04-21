# This module is used for call referenceGenomeTracker to call referenceGenomeTracker

# import pandas as pd
import numpy as np, os, csv
from ReferenceMap import refMap
from time import time
from Wig import Wig


# start = time()

# cutoffs = [50, 100, 200, 400]

# wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

# kfPath = "/home/tmhbxx3/archive/KaifuH3K4me3Fnor"


# wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"
# #
# wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]
#
# cutoffs = [x for x in range(10,310,10)]
#
# avg_size = np.zeros((len(cutoffs), len(wigFiles)))
#
# coverage = np.zeros((len(cutoffs), len(wigFiles)))
#
# peak_number = np.zeros((len(cutoffs), len(wigFiles)))
#
# for i in range(len(wigFiles)):
#     wigFile = wigFiles[i]
#     wig_obj = Wig(wigPath+wigFile)
#     for cutoff in range(10,310,10):
#         coverage[cutoff, i] = wig_obj.get_coverage(cutoff)
#         peak_number[cutoff, i] = wig_obj.get_peak_number(cutoff)
#         avg_size[cutoff, i] = coverage[cutoff, i]/peak_number[cutoff, i]
#
# np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/avg_peak_size_vs_cutoff.txt", avg_size, delimiter="\t")
# np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/coverage_vs_cutoff.txt", coverage, delimiter="\t")
# np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/peak_number_vs_cutoff.txt", peak_number, delimiter="\t")

for cutoff in [10]:
    cur_refmap = refMap(1)
    print cutoff, 'is start'
    # cur_refmap.trainMap("/home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff/pooled", cutoff=cutoff,
    #                 individual=True)

    cur_refmap.trainMap("/home/tmhbxx3/archive/WigChrSplits/code/extend_100_10/",
                        outputname='extended', cutoff=cutoff,
                        individual=False, saveRefMap=True)
    # refMap.trainMap(kfPath, surfix="regions.xls", cutoff=cutoff)
# end = time()

# print end - start

