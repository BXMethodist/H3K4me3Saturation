import numpy as np
import pandas as pd


a = pd.DataFrame([[1,2,3],[4,5,6]], columns=['A','B','C'])
print a
a['D'] = a[['A','B']].mean(axis=1)
print np.log2(4)