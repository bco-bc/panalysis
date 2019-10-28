"""
Correlation of data or vectors. May involve auto correlation. 
"""

import numpy as np
a = np.loadtxt('/localdisk/tests/dipole-moment.dat', usecols=(2,5))
print(a)
print(a.shape)
print(a.ndim)
