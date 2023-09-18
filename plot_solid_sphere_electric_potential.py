"""Plots the electric potential of a uniformly charged solid sphere, that is rho(r)=rho for r<=R and
rho(r)=0 for r>R, where R is the radius of the sphere.
"""

import logging
import sys
import numpy as np
import matplotlib.pyplot as plt

import potentials.lj_coulomb as lj_coulomb
import util


def usage():
    print()
    print('Usage: plot_solid_sphere_electric_potential.py [OPTION [value]]...')
    print('Description: Plots the electric potential due to a uniformly charged solid sphere.')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')

    dr = 0.01
    Q = 1
    R = 5.0
    r_max = 20.0
    try:
        r_max = conf['r-max']
        dr = conf['bin-size']
        R = conf['radius']
        Q = conf['charge']
    except KeyError:
        print('Using default values.')

    r = np.arange(0.0, r_max + dr, dr)
    coulomb, solid = lj_coulomb.solid_sphere_density(r, (R, Q))

    plt.plot(r/R, coulomb/solid[0], color='orange', label='Coulomb')
    plt.plot(r/R, solid/solid[0], color='blue', label='Uniformly charged solid sphere')
    plt.ylabel(r'$\psi$(r)/$\psi$(0)', style='italic')
    plt.xlabel('r/R', style='italic')
    # maximum = max(solid) + 1.0
    plt.ylim(-0.02, 2.0 * max(solid/solid[0]))
    # plt.legend()
    plt.show()
