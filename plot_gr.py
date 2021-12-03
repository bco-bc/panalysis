"""Plot radial distribution function
"""


from particle.particle_system import read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
from analysis.gr import Gr
import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    if len(sys.argv) < 5:
        print('Usage:')
        print('python fnParticleSpecs fnParticleSystem fnTrajectory spec_1 spec_2')
        raise Exception('Missing argument')

    catalog = read_particle_spec_catalog(sys.argv[1])
    particle_system = read_particle_sys(sys.argv[2], catalog)
    box = particle_system.box
    spec_1 = catalog.find(sys.argv[4])
    spec_2 = catalog.find(sys.argv[5])

    # Compute g(r)
    trajectory = Trajectory(sys.argv[3])
    analyzer = Gr(bin_size=0.01, box = box, r_max = 2.6, spec_1 = spec_1, spec_2 = spec_2)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    r, gr = analyzer.results()

    # Plot g(r)
    ones = np.ones(r.size)
    plt.plot(r, gr, color='red')
    plt.plot(r, ones, color='blue')
    plt.xlabel(r'$r$ (nm)')
    plt.ylabel(r'$g(r)$')
    plt.show()
