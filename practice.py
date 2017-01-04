import numpy as np
from scipy.stats import pearsonr
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
from time import time


a = [1,2,3]
a = np.asarray(a)

start = time()
print (a>1).sum()
end = time()

print end - start
