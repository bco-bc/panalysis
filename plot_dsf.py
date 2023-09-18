"""Plot damped shifted force potentials
"""

from potentials.gaussian_potential import GaussianPotential
from potentials.coulomb_potetial import CoulombPotential
import potentials.shifted_force as sf
import numpy as np
import matplotlib.pyplot as plt
import math

if __name__ == '__main__':
    # Assume normal distribution. 99.8% of the charge density is included at -3*sigma to 3*sigma. Assume
    # diameter of the density to be 2.0 nm = 6 * sigma (about).
    r_C = 2.0
    sigma = r_C / 6.0
    q = 1.0

    # Potentials
    gaussian_potential = GaussianPotential(sigma=sigma, q=q)
    coulomb_potential = CoulombPotential(q=q)

    # Assume nm.
    dr = 0.0001
    r = np.arange(start=dr, stop=r_C + sigma, step=dr)

    param_gaussian = {'q': -q, 'sigma': sigma}
    print(f'Cutoff distance r_C: {r_C}')
    print(f'sigma: {sigma}')
    print(f'dr: {dr}')

    # Interaction.
    r, u_r_gaussian_SF, du_dr_r_gaussian_SF = \
        sf.shifted_force(gaussian_potential, r, r_C, param=param_gaussian)
    print(f'Shifted gaussian potential at r={dr}: {u_r_gaussian_SF[0]}')
    r, u_r_coulomb_SF, du_dr_r_coulomb_SF = sf.shifted_force(coulomb_potential, r, r_C, param={'q': -q})

    # Plot potentials
    # figure = plt.figure(figsize=(8, 8))

    # Gaussian potential
    #plt.subplot(2, 1, 1)
    zeros = np.zeros(shape=len(r))
    u_0 = math.fabs(u_r_gaussian_SF[0])
    u_r_gaussian_SF /= u_0
    u_r_coulomb_SF /= u_0
    f_r_gaussian_SF = -1.0 * du_dr_r_gaussian_SF
    f_r_coulomb_SF= -1.0 * du_dr_r_coulomb_SF
    f_0 = math.fabs(f_r_gaussian_SF[0])
    f_r_gaussian_SF /= f_0
    f_r_coulomb_SF /= f_0
    plt.plot(r, u_r_gaussian_SF, color='red', label=r'$U_{G,SF}(r)$,' + rf'$\sigma={sigma}$ nm')
    plt.plot(r, u_r_coulomb_SF, color='blue', label=r'$U_{C,SF}(r)$')
    plt.plot(r, zeros, '--', color='black')
    plt.plot(r, f_r_gaussian_SF, color='green', label=r'$dU_{G,SF}(r)/dr$,' + f'$\sigma={sigma}$ nm')
    plt.plot(r, f_r_coulomb_SF, color='black', label=r'$dU_{C,SF}(r)/dr$')
    plt.xlabel(r'$r$ (nm)')
    min_du_dr = min(f_r_gaussian_SF)
    max_du_dr = max(f_r_gaussian_SF)
    min_u = min(u_r_gaussian_SF)
    max_u = max(u_r_gaussian_SF)
    #plt.ylim(1.5 * min(min_u, min_du_dr), 2.0 * max(max_u, max_du_dr))
    plt.ylim(-10.0, 2.0 * max(max_u, max_du_dr))
    plt.legend()

    # Solid sphere potential
    #plt.subplot(2, 1, 2)

    plt.show()
