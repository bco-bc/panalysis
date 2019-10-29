"""Reads time and vector from input file, computes vector's auto correlation function, and generates a plot.
"""

import sys
import numpy as np
import correlation.corr as corr
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print('Number of arguments: ', len(sys.argv))
    print('Argument List: ', str(sys.argv))
    raise Exception('Missing input filename')

fn = sys.argv[1]
print('Input file name: ', fn)

data = np.loadtxt(fn, usecols=(0, 1, 2, 3))

results = corr.acf(data)
time = results[0]
zeros = results[1]
acf = results[2]

plt.plot(time, acf)
plt.plot(time, zeros);
plt.xlabel('time (ps)')
plt.ylabel('ACF')
plt.show()
