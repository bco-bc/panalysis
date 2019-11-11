"""Plot for several approximate schemes the Coulomb interaction between two charges at given distances.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import interaction.lj_coulomb as lj_coulomb

pi = np.pi
e0 = 0.000572766
eps = 2.5           # Default value of the relative permittivity (dielectric constant)
eps_cs = eps        # Default value reaction field potential for part inside cutoff distance.
eps_rf = 78.5       # Default value reaction field potential for part outside cutoff distance.
dr = 0.1            # Distance spacing.
r0 = 0.1            # Start value distance.
alpha = 2.0         # Damping parameter, in nm^-1
kappa = 0.0         # Inverse Debye length.
q1 = -0.575         # Default charge value #1
q2 = +0.575         # Default charge value #2
rc = 2.6            # Default cutoff distance, nm.

print('Usage: python3 plot_coulomb.py (<q1> <q2> <rc> (<eps>)')
print('All arguments are optional.')
print('q1, q2 are the charge value #1, default are ', q1, ' and ', q2, ' , respectively.')
print('rc is the cutoff distance. Default is ', rc)
print('eps is relative permittivity, default value is ', eps, '.')

if len(sys.argv) == 2:
    q1 = float(sys.argv[1])
if len(sys.argv) == 3:
    q1 = float(sys.argv[1])
    q2 = float(sys.argv[2])
if len(sys.argv) == 4:
    q1 = float(sys.argv[1])
    q2 = float(sys.argv[2])
    rc = float(sys.argv[3])
if len(sys.argv) == 5:
    q1 = float(sys.argv[1])
    q2 = float(sys.argv[2])
    rc = float(sys.argv[3])
    eps = float(sys.argv[4])
    eps_cs = eps

print('q1 = ', q1)
print('q2 = ', q2)
print('rc = ', rc)
print('eps = ', eps)

C12 = 0.0
C6 = 0.0
param_exact = (q1, q2, eps, rc)
param_sf = (q1, q2, eps, rc)
param_dsf = (q1, q2, eps, rc, alpha)
param_rf = (q1, q2, eps_cs, eps_rf, rc, kappa)

r = np.arange(r0, rc, dr)
r = np.append(r, rc)
coulomb_rc = lj_coulomb.coulomb_cutoff(r, param_exact)
sf = lj_coulomb.shifted_force(r, param_sf)
dsf = lj_coulomb.damped_shifted_force(r, param_dsf)
sf3 = lj_coulomb.shifted_force_3nd_derivative(r, param_sf)
rf = lj_coulomb.reaction_field(r, param_rf)

# Plot graph.
plt.plot(r, coulomb_rc, color='black', label='RC')
label = 'DSF, alpha = ' + str(alpha) + '/nm'
plt.plot(r, dsf, color='red', label=label)
plt.plot(r, sf, color='blue', label='SF')
plt.plot(r, sf3, color='green', label='SF3')
plt.plot(r, rf, color='orange', label='RF')
s = 'q1 = ' + str(q1) + ', q2 = ' + str(q2) + ', eps = ' + str(eps) + ', rc = ' + str(rc)
plt.title(s)
plt.xlabel('r (nm)')
plt.ylabel('U(r) (kJ/mol)')
plt.legend()
plt.show()
