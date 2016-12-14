import numpy as np

a = np.array([1,1,-1,-2,-3,4,5])
asign = np.sign(a)


signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
print signchange

itemindex = np.where(signchange==1)

print itemindex