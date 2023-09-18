"""Plot distance between two particles
"""

import logging
import sys

import matplotlib.pyplot as plt
import numpy as np

import util
from analysis import analyze
from analysis.distance import Distance
from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
from simulation.trajectory import Trajectory


def usage():
    print()
    print('Usage: python plot_distance.py [OPTION [value]]...')
    print('Description: Computes and plots the distance between two specific particles, '
          'identified by their particle identifiers, or between particles of given specifications in groups. '
          'Either ids or specs must be provided, not both.')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')

    print()
    print('Optional arguments:')
    print('-id-1 ID: Particle identifier particle #1.')
    print('-id-2 ID: Particle identifier particle #2.')
    print('-spec-1 ID: Particle specification particle #1.')
    print('-spec-2 ID: Particle specification particle #2.')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print('-dt V: V is the value of the time difference between entries in the trajectory. Default is 0.002.')
    print('-pbc-1 V: V is the direction along which PBC should be applied. One of {x, y, z}. Default is to apply '
          'PBC in all directions.')
    print('-n-skip VALUE: Skip the first VALUE state(s) from trajectory. Default is 0.')
    print('-R V: V is a reference distance. Default is 0')
    print()


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        id_1 = conf['id-1']
        id_2 = conf['id-2']
        spec_1 = conf['spec-1']
        spec_2 = conf['spec-2']
        fn_trajectory = conf['fn-trajectory']
        pbc = conf['pbc']
        dt = conf['dt']
        n_skip = int(conf['n-skip'])
        distance = float(conf['distance'])
    except (KeyError, FileNotFoundError, FileExistsError):
        usage()
        sys.exit("Missing input argument(s) and/or non-existing files.")

    # Monitor distance
    specs = len(spec_1) and len(spec_2)
    print(f'PBC: {pbc}')
    trajectory = Trajectory(fn_trajectory)
    if not specs:
        analyzer = Distance(dt=dt, id_1=id_1, id_2=id_2, pbc=pbc)
    else:
        analyzer = Distance(dt=dt, spec_1=spec_1, spec_2=spec_2, pbc=pbc)
    analyze.perform(analyzer, particle_system, trajectory, n_skip)
    trajectory.close()
    if not specs:
        t, dis, x, y, z = analyzer.results()
        R = np.zeros(shape=len(t))
        for k in np.arange(0, len(t)):
            R[k] = distance
    else:
        dis, x, y, z = analyzer.results()
        R = np.zeros(shape=len(dis))
    if len(dis) == 0:
        print('No particle pairs found.')
        quit(0)

    print(f'{len(dis)}: Length of dis')
    ave_dis = sum(dis[:]) / len(dis)
    print(f'Average distance: {ave_dis}')

    # Plot results
    figure = plt.figure()

    # Plot distance versus time
    plt.subplot(3, 1, 1)
    if not specs:
        t = np.arange(0, len(dis))
        plt.plot(t, dis, color='blue', label='r')
        plt.plot(t, x, color='grey', linestyle='dashdot', label='x')
        plt.plot(t, y, color='grey', linestyle='--', label='y')
        plt.plot(t, z, color='grey', linestyle='solid', label='z')
        plt.plot(t, R, color='red', label=r'R$_{ref}$')
        plt.xlabel(r't')
        plt.ylabel(r'$r(t)$')
        plt.legend()

    # Plot coordinate densities
    plt.subplot(3, 1, 2)
    hist_x, bin_edges_x = np.histogram(x, density=True)
    n = len(bin_edges_x) - 1
    plt.plot(bin_edges_x[0:n], hist_x, color='red', label='x')
    hist_y, bin_edges_y = np.histogram(y, density=True)
    plt.plot(bin_edges_y[0:n], hist_y, color='blue', label='y')
    hist_z, bin_edges_z = np.histogram(z, density=True)
    plt.plot(bin_edges_z[0:n], hist_z, color='green', label='z')
    plt.xlabel(r'x')
    plt.ylabel(r'p(x)')
    plt.legend()

    plt.subplot(3, 1, 3)
    # p_dis, p_dis_edges = np.histogram(a=dis, range=(0.0, 1.0), density=True)
    p_dis, p_dis_edges = np.histogram(a=dis, density=True)
    m = len(p_dis_edges) - 1
    plt.plot(p_dis_edges[0:n], p_dis, color='red', label='p(r)')
    plt.xlabel(r'r')
    plt.ylabel(r'p(r)')

    # Probability for being an interval greater than threshold value (reference distance)
    probability = 0
    for k in np.arange(0, len(p_dis)-1):
        if (p_dis_edges[k]) > distance:
            probability += p_dis[k] * (p_dis_edges[k+1] - p_dis_edges[k])
    print(f'Probability(r > reference distance): {probability}')
    figure.tight_layout(pad=1.0)
    plt.show()
