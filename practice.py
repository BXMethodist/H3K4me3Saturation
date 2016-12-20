import numpy as np
import pandas as pd

a = np.array([1,1,-1,-2,-3,4,5, 6])
# asign = np.sign(a)
#
# a = a.reshape((2,4))
#
# print a, type(a)
#
# signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
# print signchange
#
# itemindex = np.where(signchange==1)[0]

# print type(itemindex), itemindex

a[a<0] = 0
print a
