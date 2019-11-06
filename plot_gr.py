import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print('Number of arguments: ', len(sys.argv))
    print('Argument List: ', str(sys.argv))
    print('Usage: python3 plot_gr <filename>')
    print(' filename: Data input file. Must be provided.')
    raise Exception('Missing input filename')

fn = sys.argv[1]
print('Input file name: ', fn)

# Read data.
data = np.loadtxt(fn, usecols=(0, 1))

r = data[:, 0]  # Each row, first column.
gr = data[:, 1] # Each rown, second column.
ones = np.ones((r.size))

plt.plot(r, gr)
plt.plot(r, ones)
plt.xlabel('r (nm)')
plt.ylabel('g(r)')
plt.legend()
plt.show()
