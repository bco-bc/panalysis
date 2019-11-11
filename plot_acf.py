"""Reads time and values of scalar quantity (components of vector) from input file, computes quantity' (vector's)
auto correlation function (ACF), and generates a plot. Input file should contain 2 or more columns of data. The first
column represents time, the remaining one(s) holds (hold) the values of the scalar quantity (the components
of the vector).
"""

import sys
import numpy as np
import correlation.corr as corr
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print('Number of arguments: ', len(sys.argv))
    print('Argument List: ', str(sys.argv))
    print('Usage: python3 v_acf.py <filename> <n_use>')
    print(' filename: Data input file. Must be provided.')
    print(' n_use: Number of elements to use from ACF. Default is 50.')
    raise Exception('Missing input filename')

fn = sys.argv[1]
print('Input file name: ', fn)

n_use = 50;
if len(sys.argv) == 3:
    n_use = int(sys.argv[2])

print('Using', n_use, 'elements of ACF.')

# Read data
data = np.loadtxt(fn, usecols=(0, 1, 2, 3))

# Extract time and vector.
results = corr.acf(data, n_use)
time = results[0]
zeros = np.zeros(n_use)
acf = results[1]

plt.plot(time, acf)
plt.plot(time, zeros);
plt.xlabel('time (ps)')
plt.ylabel('ACF')
plt.show()
