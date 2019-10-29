"""
Correlation of data or vectors.
"""

import numpy as np
import tidynamics as td

a = np.loadtxt('/localdisk/tests/dipole-moment.dat', usecols=(0,1,2,3))
print(a)
print(a.shape)
print(a.ndim)
