"""Plots the reaction potential inside a sphere
"""

import logging
import matplotlib.pyplot as plt
import numpy as np

import polarization.reaction_potential as reaction_potential


def usage():
    print()
    print('Usage: python plot_reaction_potential')


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    param = {'radius': 3.0, 'epsilon': 78.5, 'bin-size': 0.02}

    bin_size = param['bin-size']
    radius = param['radius']
    n = int(radius / bin_size)
    bin_size = radius / n
    Q = 1.0
    param['bin-size'] = bin_size

    logging.info('Parameters: ' + str(param))
    z, r_p = reaction_potential.image_approximation_q(Q, param)

    fn = 'reaction-potentials-im.dat'
    f = open(fn, 'w', encoding='utf-8')
    for i in np.arange(0, len(z)):
        z_i = str(z[i])
        r_p_i = str(r_p[i])
        s = z_i + ' ' + r_p_i + '\n'
        f.writelines(s)
    f.close()

    plt.plot(z, r_p, color='red')
    plt.xlabel(r'$r$ (nm)')
    plt.ylabel(r'$\phi_R$ (r)')
    plt.show()
