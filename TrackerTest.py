# import pandas as pd
import numpy as np, os, csv
from referenceGenomeTracker import *
from time import time

refMap = H3K4me3Saturation(5)
# start = time()
refMap.trainMap("./KaifuH3K4me3Fnor", cutoff=300)
# end = time()

# print end - start

