"""Plots the total electric current autocorrelation function over given time interval.
"""

from particle.particle_system import read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import ParticleSystem
from simulation.trajectory import Trajectory
from analysis.current_acf import CurrentACF
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
import logging
import util


def usage():
    print()
    print('Usage: python plot_current_acf.py [OPTION [value]]...')
    print('Description: Plots total electric current autocorrelation function.')
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
    print('-E0 VALUE: VALUE is the strength of the external electric field (in V/nm).')
    print('-direction-EO VALUE: VALUE is direction of the external field, one of \'x\', \'y\', or \'z\'.')
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
    ef_0 = np.zeros(shape=3)
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        dt = conf['dt']
        t_max = conf['t-max']
        temperature = conf['temperature']
        fn_trajectory = conf['fn-trajectory']
        e_0_strength = conf['E0']
        d = conf['direction-E0']
        if d == 'x':
            ef_0[0] = e_0_strength
        elif d == 'y':
            ef_0[1] = e_0_strength
        else:
            ef_0[2] = e_0_strength
    except (KeyError, FileExistsError, FileNotFoundError):
        usage()
        sys.exit("Missing input argument(s).")

    # Compute current ACF.
    trajectory = Trajectory(fn_trajectory)
    analyzer = CurrentACF(dt=dt, t_max=t_max)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    t, current_acf = analyzer.results()
    logging.info(f'Length of time interval: {analyzer.t_max}')
    logging.info(f'Number of states/entries in trajectory: {analyzer.counter}')
    logging.info(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')

    # Cubic interpolation.
    if np.linalg.norm(ef_0) > 0.0:
        f_current_acf = interpolate.interp1d(t, current_acf, kind='cubic')

        # Function for integration.
        def current_acf_f(s: float) -> float:
            return f_current_acf(s)

        t_new = np.linspace(0, t[t.size-1], num=2*t.size)
        ave_current = np.zeros(shape=(t_new.size, 3))
        t_index = 0
        val, err = integrate.quadrature(func=current_acf_f, a=0.0, b=t_new[t_new.size-1], maxiter=1000)
        logging.info(f'Value of integration of <J(t)J(0)>: {val}. Error= {err}')
        for time in t_new:
            ave_current[t_index] = ef_0 * val
            t_index += 1
    else:
        t_new = np.zeros(shape=0)
        ave_current = np.zeros(shape=0)

    # Plot electric current ACF
    zeros = np.zeros(shape=analyzer.n_bins)
    plt.plot(t, current_acf, color='red')
    plt.plot(t, zeros, '--', color='black')
    plt.xlabel(r'$t$ (ps)')
    plt.ylabel(r'$\frac{<J(t)J(0)>}{<J(0)J(0)>}$ ((e.nm/ps)$^2$)')

    if np.linalg.norm(ef_0) > 0.0:
        ave = np.zeros(shape=t_new.size)
        a_index = 0
        for a in ave_current:
            ave[a_index] = np.linalg.norm(a)
            a_index += 1
        plt.plot(t_new, ave, color='blue')

    plt.show()
