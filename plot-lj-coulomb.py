import sys
import numpy as np
import matplotlib.pyplot as plt
import interaction.lj_coulomb as lj_coulomb

pi = np.pi
e0 = 0.000572766
eps = 2.5           # Default value of the relative permittivity (dielectric constant)
dr = 0.01           # Distance spacing.

if len(sys.argv) < 5:
    print('Number of arguments: ', len(sys.argv))
    print('Argument List: ', str(sys.argv))
    print('Usage: python3 plot_lj_coulomb <q1> <q2> <C12> <C6> (<eps>)')
    print('Use \'molecular units\'. Relative permittivity \'eps\' is optional, default value is 2.5.')
    raise Exception("Missing interaction parameters.")

if len(sys.argv) == 6:
    eps = float(sys.argv[5])

q1 = float(sys.argv[1])
q2 = float(sys.argv[2])
C12 = float(sys.argv[3])
C6 = float(sys.argv[4])
param = (q1, q2, eps, C12, C6)

# Calculate interaction.
r = np.arange(0.4, 4.0, dr)
results = lj_coulomb.potential(r, param)
total = results[0]
el = results[1]
lj = results[2]

# Plot graph.
plt.plot(r, total, color='black', label='total')
plt.plot(r, lj, color='red', label='lj')
plt.plot(r, el, color='blue', label='el')
plt.xlabel('r (nm)')
plt.ylabel('U(r)')
plt.show()
