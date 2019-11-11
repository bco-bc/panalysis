"""Extracts the temperature from a simulation data output file, and generates a plot versus time, along with
requested reference temperature and the time average.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print('Number of arguments: ', len(sys.argv))
    print('Argument List: ', str(sys.argv))
    print('Usage: python3 plot_gr.py <filename> (<Tref>)')
    print(' filename: Data input file. Must be provided.')
    print(' Tref : Reference temperature. If ommited, 298.15 K is assumed.')
    raise Exception('Missing input filename')

fn = sys.argv[1]
print('Input file name: ', fn)
Tref = 298.15
if len(sys.argv) == 3:
    Tref = float(sys.argv[2])

print('Reference temperature: ', Tref)

# Read data.
data = np.loadtxt(fn, usecols=(1, 6))
t = data[:, 0]
T = data[:, 1]
T_ref = np.ones(t.size) * Tref
ave = np.arange(0)
counter = 0
sum = 0.0
for temp in T:
    counter += 1
    sum += temp
    v = sum / counter
    ave = np.append(ave, v)


# Plot
a = ave[ave.size - 1]
plt.plot(t,T, label='T')
label = '<T> = ' + str(a)
plt.plot(t, ave, label=label)
label = 'Tref = ' + str(Tref)
plt.plot(t, T_ref, label=label)
plt.xlabel('t (ps)')
plt.ylabel('T (K)')
plt.legend()
plt.show()


