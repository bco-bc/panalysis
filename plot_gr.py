"""Plot radial distribution function
"""


from particle.particle_system import ParticleSystem, read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
from analysis.gr import Gr
import sys
import numpy as np
import matplotlib.pyplot as plt
import util
import logging


def usage():
    print()
    print('Usage: python plot_gr.py [OPTION [value]]...')
    print('Description: Computes and plots the radial distribution function g(r) over '
          'given distance interval.')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('-r-max VALUE: VALUE is the length of distance interval (in nm) over which g(r) is calculated.')
    print('-spec-1 NAME: NAME is one of two particle specification names.')
    print('-spec-2 NAME: NAME is one of two particle specification names.')
    print()
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print('-pbc-1 DIRECTION: DIRECTION is the direction along PBD should be applied.')
    print()


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')
    fn_trajectory = 'trajectory.dat'
    spec_1 = None
    spec_2 = None
    r_max = 2.6
    particle_system = ParticleSystem()
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        spec_1 = catalog.find(conf['spec-1'])
        spec_2 = catalog.find(conf['spec-2'])
        fn_trajectory = conf['fn-trajectory']
        r_max = conf['r-max']
        pbc = conf['pbc']
    except (KeyError, FileNotFoundError, FileExistsError):
        usage()
        sys.exit("Missing input argument(s).")

    # Compute g(r)
    print(f'PBC: {pbc}')
    trajectory = Trajectory(fn_trajectory)
    analyzer = Gr(bin_size=0.01, box=particle_system.box, r_max=r_max, spec_1=spec_1, spec_2=spec_2, pbc = pbc)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    r, gr = analyzer.results()

    # Plot g(r)
    ones = np.ones(r.size)
    plt.plot(r, gr, '-', color='red')
    plt.plot(r, ones, '--', color='blue')
    plt.xlabel(r'$r$ (nm)')
    plt.ylabel(r'$g(r)$')
    plt.show()
