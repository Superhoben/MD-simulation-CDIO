import numpy as np
from scipy.spatial.distance import cdist

a = np.array([[1,2,3], [1,2,4], [2,3,4]])
b = cdist(a,a)

c = np.delete(b[0], 0)
cur_min = min(c)

