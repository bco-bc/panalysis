"""Plot number density in slices
"""

from particle.particle_system import read_particle_sys, ParticleSystem
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
from analysis.slice_number_density import SliceNumberDensity
import matplotlib.pyplot as plt
import numpy as np
import sys
import util
import logging

def usage():
    print()
    print('Usage: python plot_gr.py [OPTION [value]]...')
    print('Description: Computes and plots the particle number density in slices.')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('-spec NAME: NAME is the particle specification name.')
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')
    fn_trajectory = 'trajectory.dat'
    spec = None
    particle_system = ParticleSystem()
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        spec = catalog.find(conf['spec'])
    except (KeyError, FileNotFoundError, FileExistsError):
        usage()
        sys.exit("Missing input argument(s).")

    # Compute number density in slices
    trajectory = Trajectory(fn_trajectory)
    analyzer = SliceNumberDensity(0.01, particle_system.box, spec)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    z, number_densities = analyzer.results()

    n_particles = 0
    for p in particle_system.all:
        if p.spec.name is spec.name:
            n_particles += 1
    volume = particle_system.box_volume()
    bulk_number_densities = np.zeros(shape = len(number_densities))
    k_values = np.arange(0, len(number_densities))
    for k in k_values:
        bulk_number_densities[k] = n_particles / volume

    # Plot
    plt.plot(z, number_densities, color='red')
    plt.plot(z, bulk_number_densities, color='blue')
    plt.xlabel(r'$z$ (nm)')
    plt.ylabel(r'$\rho(z)$')
    plt.show()
