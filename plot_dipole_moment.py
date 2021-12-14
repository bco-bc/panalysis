"""Plots total dipole moment autocorrelation function over given time interval.
"""

import logging
import sys
import util
from simulation.trajectory import Trajectory
from particle.particle_system import read_particle_sys, ParticleSystem
from particle.particle_spec_catalog import read_particle_spec_catalog
from analysis.dipole_moment import DipoleMoment
import simulation.trajectory
import matplotlib.pyplot as plt
import numpy as np

def usage():
    print()
    print('Usage: python plot_dipole_moment.py [OPTION [value]]...')
    print('Description: Plots total dipole moment autocorrelation function.')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('-dt VALUE: VALUE is the time difference (in ps) between states in trajectory.')
    print('-t-max VALUE: VALUE is the length of time interval (in ps) over which VAF is calculated.')
    print()
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print('-T VALUE: VALUE is the temperature (in K). Default is 298.15 K.')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)
    conf = util.parse_argv.parse(sys.argv)
    logging.info(f'Received program arguments {conf}')
    fn_trajectory = 'trajectory.dat'
    temperature = 298.15
    t_max = 0.0
    spec = ''
    dt = 0.0
    particle_system = ParticleSystem()
    cos_a_bin_size = 0.01
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        dt = conf['dt']
        t_max = conf['t-max']
        temperature = conf['temperature']
        fn_trajectory = conf['fn-trajectory']
    except (KeyError, FileExistsError, FileNotFoundError):
        usage()
        sys.exit("Missing input argument(s).")

    # Compute dipole moment.
    trajectory = Trajectory(fn_trajectory)
    analyzer = DipoleMoment(dt=dt, t_max=t_max, cos_a_bin_size=cos_a_bin_size)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    t, m, ave_m, p = analyzer.results()
    print(f'Dipole moment averaged of trajectory: {ave_m}')
    print(f'Norm average dipole moment: {np.linalg.norm(ave_m)}')

    # Plot results.
    figure = plt.figure(figsize=(8, 8))

    plt.subplot(2, 1, 1)
    t /= 1000.0    # To nanoseconds.
    norm_m = np.zeros(shape=len(m))
    norm_m_t = np.zeros(shape=len(m))
    m_t = np.array([0.0, 0.0, 0.0])
    for i in np.arange(0, len(m)):
        m_t += m[i]
        norm_m[i] = np.linalg.norm(m[i])
        ave_m_t = m_t / (1.0 + i)
        norm_m_t[i] = np.linalg.norm(ave_m_t)
    zeros = np.zeros(shape=t.shape)
    plt.plot(t, m, color='grey')
    plt.plot(t, norm_m_t, color='red')
    plt.plot(t, zeros, '--', color='black')
    plt.xlabel(r'$t$ (ns)')
    plt.ylabel(r'$M(t)_k, k \in{x,y,z}$')

    plt.subplot(2, 1, 2)
    cos_a = np.arange(-1.0, +1.0, step=analyzer.cos_a_bin_size)
    plt.plot(cos_a, p, color='blue')
    plt.xlabel(r'$\cos(\alpha)$')
    plt.ylabel(r'$p(\cos(\alpha))$')
    plt.ylim(0.0, 1.0)

    figure.tight_layout(pad=1.0)
    plt.show()
