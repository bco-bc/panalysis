"""Plot velocity profile in a channel in given direction.
"""

import logging
import sys

import matplotlib.pyplot as plt
import numpy as np
import util
from analysis.velocity_profile import VelocityProfile
from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys, ParticleSystem
from simulation.trajectory import Trajectory


def usage():
    print()
    print('Usage: python plot_velocity_profile [OPTION [value]]...')
    print("Computes and plots velocity profile in slice in a given direction.")
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print('-d DIRECTION: DIRECTION is the direction along which the velocity profile is computed. '
          'Default is the z-direction.')
    print('-d2 DIRECTION: DIRECTION is the second or other direction. Default is the x-direction.')
    print('-bin-size VALUE: VALUE is the bin size. Default is 0.1')
    print('-pbc-1 V: V is the direction along which PBC should be applied. One of {x, y, z}. Default is to apply '
          'PBC in all directions.')
    print('-exclude SPEC: SPEC is the particle specification names of particles that must be '
          'excluded. This may be one name or a comma-separated list of names (no spaces)')


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')
    particle_system = ParticleSystem()
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        bin_size = float(conf['bin-size'])
        direction = conf['direction']
        direction_2 = conf['second-direction']
        pbc = conf['pbc']
        exclude = conf['exclude']
    except (KeyError, FileNotFoundError, FileExistsError):
        usage()
        sys.exit("Missing input argument(s).")

    print(f'Calculating velocity profile in the {direction}-direction')
    trajectory = Trajectory(conf['fn-trajectory'])
    analyzer = VelocityProfile(bin_size=bin_size,
                               box=particle_system.box,
                               direction=direction,
                               direction_2=direction_2,
                               pbc=pbc,
                               exclude=exclude)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    r, velocity_profile = analyzer.results()
    zeros = np.zeros(shape=len(r))
    plt.plot(r, velocity_profile, color='red')
    plt.plot(r, zeros, '--', color='black')
    plt.xlabel(direction_2)
    plt.ylabel(f'$<v_{direction}>$')
    plt.show()
