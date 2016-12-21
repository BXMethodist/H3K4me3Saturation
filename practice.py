import numpy as np
from scipy.stats import pearsonr
import pandas as pd

a = [0,1,1, 0]
b = [1,2,2,1]
c =[ -1, -2, -2, -1]
d = [ 0, 0, 1, 1]


result = []
result.append(a)
result.append(b)
result.append(c)
result.append(d)


result = np.asarray(result, dtype=np.float64)

r = result.copy()
r[0,0] = np.nan

print np.mean(result, axis = 0)

# print result
#
# print r


pe = np.corrcoef(result)
#
pe[0,0] = np.nan

index = np.argmin(pe)

print index

print pe

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
