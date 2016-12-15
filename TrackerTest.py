# This module is used for call referenceGenomeTracker to call referenceGenomeTracker

# import pandas as pd
import numpy as np, os, csv
from refMap import *
from time import time

refMap = refMap(5)
# start = time()

cutoffs = [25,50]
for cutoff in cutoffs:
    refMap.trainMap("/home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff/pooled", cutoff=cutoff)
# end = time()

# print end - start

