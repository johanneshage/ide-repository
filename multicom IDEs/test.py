from collections import defaultdict
import scipy.sparse
import functools
import numpy as np


I = 2
m= 20

# def_dict = np.dtype(functools.partial(defaultdict, list))
fpi = defaultdict(list)
# fp = lil_matrix((I, m), dtype=functools.partial(defaultdict, list))

fpi[-1]

eind = 2
print("defdict0", fpi[eind])
fpi[eind].append(5)
print("defdict", fpi[eind])

#fp[0,0] = list([])
#fp[0, 0].append((0, 1))
#fp[1, 11].append((1, 2))
#fp[0,0].append((1.1, 0))

#print("00", fp[0,0])
#fp[0 , 0] += [1]
#np.append(fp[0, 0], 1)

#print("FP", fp)

liste1 = [0, 1.5, 0, 0, 1.9, 0]
matrix1 = scipy.sparse.csr_matrix(liste1)
print(matrix1)

for e in range(6):
    print("LISTE", matrix1[0, e])

bounds = set()
print("bounds1", bounds)
bounds += {1}
print("bounds2", bounds)
bounds += {2}
print("bounds3", bounds)