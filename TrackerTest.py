# import pandas as pd
import numpy as np, os, csv
from referenceGenomeTracker import *
from time import time

refMap = H3K4me3Saturation()
# start = time()
refMap.trainMap("./KaifuH3K4me3Fnor")
# end = time()

# print end - start

# print refMap.results
