"""Plots the electric potential of a Gaussian charge density given by rho(r)=A*exp(-ar^2), where A and a
are constants.
"""

import logging
import sys
import numpy as np
import matplotlib.pyplot as plt

import potentials.lj_coulomb as lj_coulomb
import util


def usage():
    print()
    print('Usage: plot_gaussian_electric_potential.py [OPTION [value]]...')
    print('Description: Plots the electric potential due to a Gaussian charge density.')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')
    r_max = 30.0
    dr = 0.001
    sigma = 1.0
    Q = 1
    try:
        r_max = conf['r-max']
        dr = conf['bin-size']
    except KeyError:
        print('Using default values.')

    r = np.arange(0.0, r_max + dr, dr)
    coulomb, gaussian = lj_coulomb.gaussian_density(r, (sigma, Q))

    plt.plot(r/sigma, coulomb/gaussian[0], color='orange', label='Coulomb')
    plt.plot(r/sigma, gaussian/gaussian[0], color='blue', label='Gaussian')
    plt.ylabel(r'$\psi$(r)/$\psi$(0)', style='italic')
    plt.xlabel(r'r/$\sigma$', style='italic')
    plt.ylim(-0.02, 2.0 * max(gaussian/gaussian[0]))
    plt.legend()
    plt.show()
