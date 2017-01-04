# This module is used for call referenceGenomeTracker to call referenceGenomeTracker

# import pandas as pd
import numpy as np, os, csv
from refMap import *
from time import time

refMap = refMap(1)
# start = time()

cutoffs = [50, 100, 200, 400]

# wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

kfPath = "/home/tmhbxx3/archive/KaifuH3K4me3Fnor"

for cutoff in cutoffs:
    refMap.trainMap("/home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff/pooled", cutoff=cutoff)
    # refMap.trainMap(kfPath, surfix="regions.xls", cutoff=cutoff)
# end = time()

# print end - start

