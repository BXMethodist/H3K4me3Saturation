import numpy as np

a = np.asarray([0,1,1,1,0,1,1,1,0])

b = np.where(((np.roll(a, 1) - a) != 0).astype(int)==1)

print b
c = b - np.roll(b, 1)

print a[1::2]

a = np.insert(a, 5, 0)

print a