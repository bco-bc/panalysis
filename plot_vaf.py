"""Plot velocity autocorrelation function (VAF) over a given time interval.
"""

from analysis.vaf import VAF
from particle.particle_system import read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
import scipy.constants as constants
import logging
import util


def usage():
    print()
    print('Usage: python plot_vaf.py [OPTION [value]]...')
    print('Description: Computes and plots the velocity autocorrelation function (VAF).')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('-dt VALUE: VALUE is the time difference between states in trajectory.')
    print('-t-max VALUE: VALUE is the length of time interval (in ps) over which VAF is calculated.')
    print('-spec NAME: NAME is the particle specification name.')
    print()
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print('-T VALUE: VALUE is the temperature (in K). Default is 298.15 K.')
    print()


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    fn_trajectory = 'trajectory.dat'
    temperature = 298.15
    t_max = 0.0
    spec = ''
    dt = 0.0
    particle_system = ''
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        spec = catalog.find(conf['spec'])
        dt = conf['dt']
        t_max = conf['t-max']
        temperature = conf['temperature']
        fn_trajectory = conf['fn-trajectory']
    except (KeyError, FileNotFoundError, FileExistsError):
        print(f'Received arguments: {conf}')
        usage()
        sys.exit("Missing input argument(s).")

    # Compute VAF.
    trajectory = Trajectory(fn_trajectory)
    bin_size = dt
    analyzer = VAF(dt, t_max, spec)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    t, vaf = analyzer.results()

    logging.info(f'Number of particles with particle specification \'{analyzer.spec.name}\': {analyzer.n_specs}')
    logging.info(f'Length of time interval: {analyzer.t_max}')
    logging.info(f'Number of states/entries in trajectory: {analyzer.counter}')
    logging.info(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')

    # Cubic interpolation.
    f_vaf = interpolate.interp1d(t, vaf, kind='cubic')
    n = t.size * 2
    t_new = np.linspace(0, t[t.size - 1], num=2*t.size, endpoint=True)
    vaf_new = f_vaf(t_new)

    # Function for integration.
    def vaf_f(time) -> float:
        return f_vaf(time)

    a = 0.0
    b = t[t.size-1]
    val, err = integrate.quadrature(func=vaf_f, a=a, b=b, maxiter=1000)
    k_B = constants.k * constants.N_A / 1.0e+03  # Boltzmann constant In kJ/(mol K)
    D = k_B * temperature / spec.mass * val
    print(f'Diffusion (self) constant: {D} nm^2/ps')

    # Plot VAF
    zeros = np.zeros(shape=analyzer.n_bins)
    plt.plot(t, vaf, 'o', color='red')
    plt.plot(t_new, vaf_new, '--', color='blue')
    plt.xlabel(r'$t$ (ps)')
    plt.plot(t, zeros, '--', color='black')
    plt.ylabel(r'$<v(t+s)v(t)$')
    plt.legend(['data', 'cubic'])
    plt.show()
