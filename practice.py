import numpy as np

a = [ 1, 2, 3]
a = np.sign(a)
print ((np.roll(a, 1) - a) != 0).astype(int)