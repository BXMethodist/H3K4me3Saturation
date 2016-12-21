import numpy as np
from scipy.stats import pearsonr

a = [0,1,1, 0]
b = [1,2,2,1]
c =[ -1, -2, -2, -1]
d = [ 0, 0, 1, 1]


result = []
result.append(a)
result.append(b)
result.append(c)
result.append(d)


result = np.asarray(result)


pe = np.corrcoef(result)

index = np.argmin(pe)

print index

min = pe.flat[index]

print pe, min

print np.where(pe==min)
