import numpy as np
from scipy.stats import pearsonr
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd

a = np.asarray([1,2,3])
print a.shape[0]


# a = np.asarray([[1,2,3,4,5,6,7,8,9,10,11,12], [1,2,3,4,5,6,7,8,9,10,11,12], [1,2,3,4,5,6,7,8,9,10,11,12]])
#
# a = a *[1,2,3]
#
# print a

# df = pd.DataFrame([[1,2,3,4,5,6,7,8,9,10,11,12], [1,2,3,4,5,6,7,8,9,10,11,12], [1,2,3,4,5,6,7,8,9,10,11,12]])
#
# # b = pd.rolling_mean(df, 3, center=True, min_periods=1)
#
# print df.rolling(window=2, center=True, axis=1).mean().values
#
# x1 = np.arange(12).reshape((3, 4))
# x2 = np.asarray([1,2,3]).reshape(3,1)
# print x1.shape
#
# print x2, x2.shape
#
# print np.multiply(x1, x2)


# a = [1,1,0, 0]
# b = [1,2,2,1][1,2,3,4,5,6,7,8,9,10,11,12]
# c =[ -1,-2,-2,-1]
# d = [ 0, 0, 1, 1]
#
#
# result = []
# result.append(a)
# result.append(b)
# result.append(c)
# result.append(d)
#
# print np.argmax(result, axis=0)
#
#
# result = np.asarray(result, dtype=np.float64)
#
# print cosine_similarity(result)

# pe = np.corrcoef(result)
#
# print pe
# #
# pe[0,0] = np.nan
#
# index = np.argmin(pe)
#
# print index
#
# print pe

# index = np.argmin(pe)
#
# print index
#
# min = pe.flat[index]
#
# print pe, min
#
# r = np.asarray(np.where(pe==min))
#
# r = r.T
#
#
# distinct_set = set()
#
# for i in range(r.shape[0]):
#     x, y = r[i]
#     if (x, y) not in distinct_set and (y, x) not in distinct_set:
#         distinct_set.add((x, y))
# distinct_set = list(distinct_set)
# print distinct_set
#
# data = pd.DataFrame(result)
# print data
#
# print data.ix[0,:].values
#
# print len(np.zeros(10))
