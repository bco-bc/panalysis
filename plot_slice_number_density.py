"""Plot number density in slices
"""

import logging
import sys

import matplotlib.pyplot as plt

import util
from analysis.slice_number_density import SliceNumberDensity
from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys, ParticleSystem
from simulation.trajectory import Trajectory


def usage():
    print()
    print('Usage: python plot_slice_number_density.py [OPTION [value]]...')
    print('Description: Computes and plots the particle number density in slices in a given direction.')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('-spec NAME: NAME is the particle specification name.')
    print('-r-max VALUE: VALUE is the length of the distance interval for slices.')
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print('-d DIRECTION: DIRECTION is the direction along which the number density in slices is computed. '
          'Default is the z-direction.')
    print('-bin-size VALUE: VALUE is the bin size for the particle number density.')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')
    fn_trajectory = 'trajectory.dat'
    spec = None
    direction = 'z'
    particle_system = ParticleSystem()
    r_max = 0.0
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        spec = catalog.find(conf['spec'])
        direction = conf['direction']
        r_max = conf['r-max']
        pbc = conf['pbc']
        bin_size = float(conf['bin-size'])
    except (KeyError, FileNotFoundError, FileExistsError):
        usage()
        sys.exit("Missing input argument(s).")

    # Compute number density in slices
    print(f'Calculating probability density function for slices in the {direction} direction for {spec.name}.')
    trajectory = Trajectory(fn_trajectory)
    analyzer = SliceNumberDensity(bin_size=bin_size, box=particle_system.box, r_max=r_max, spec=spec, direction=direction, pbc=pbc)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    r, pdf = analyzer.results()

    # Plot
    plt.plot(r, pdf, color='red')
    plt.xlabel(r'$r$ (nm)')
    plt.ylabel(r'$p(r)$')
    plt.show()
