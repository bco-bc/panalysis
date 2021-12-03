"""Plot number density in slices
"""

from particle.particle_system import read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
from analysis.slice_number_density import SliceNumberDensity
import matplotlib.pyplot as plt
import numpy as np
import sys

if __name__ == '__main__':

    if len(sys.argv) < 4:
        print('Usage:')
        print('python fnParticleSpecs fnParticleSystem fnTrajectory spec')
        raise ValueError('Missing arguments')

    catalog = read_particle_spec_catalog(sys.argv[1])
    particle_system = read_particle_sys(sys.argv[2], catalog)
    box = particle_system.box
    spec = catalog.find(sys.argv[4])

    # Compute number density in slices
    trajectory = Trajectory(sys.argv[3])
    analyzer = SliceNumberDensity(0.01, box, spec)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    z, number_densities = analyzer.results()

    n_particles = 0
    for p in particle_system.all:
        if p.spec.name is spec.name:
            n_particles += 1
    volume = box[0] * box[1] * box[2]
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
