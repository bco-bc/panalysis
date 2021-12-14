import numpy as np

from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory
from analysis.dipole_moment import DipoleMoment
import matplotlib.pyplot as plt
import logging

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)

    fn = '/home/juffer/simulations/electrolyte/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)
    fn = '/home/juffer/simulations/electrolyte/trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)

    analyzer = DipoleMoment(dt=20.0, t_max=500, cos_a_bin_size=0.01)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    t, m, ave_m, p = analyzer.results()

    # To nanoseconds.
    t /= 1000

    print(f'Dipole moment averaged of trajectory: {ave_m}')
    print(f'Norm average dipole moment: {np.linalg.norm(ave_m)}')
    print(f'Sum of probabilities, P(cos(a)) * d(cos(a)): {sum(p[:]) * analyzer.cos_a_bin_size}')

    norm_m = np.zeros(shape=len(m))
    for i in np.arange(0, len(m)):
        norm_m[i] = np.linalg.norm(m[i])

    figure = plt.figure(figsize=(8, 8))
    plt.subplot(2, 1, 1)
    zeros = np.zeros(shape=t.shape)
    plt.plot(t, m, color='grey')
    plt.plot(t, zeros, '--', color='black')
    plt.xlabel(r'$t$ (ns)')
    # plt.plot(t, norm_m, color='red')
    plt.subplot(2, 1, 2)
    cos_a = np.arange(-1.0, +1.0, step=analyzer.cos_a_bin_size)
    plt.plot(cos_a, p, color='blue')
    plt.ylim(0.0, 1.0)
    plt.show()
