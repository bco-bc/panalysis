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
    t, current, current_acf = analyzer.results()
    logging.info(f'Length of time interval: {analyzer.t_max}')
    logging.info(f'Number of states/entries in trajectory: {analyzer.counter}')
    logging.info(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')
    logging.info(f'Length of trajectory: {analyzer.counter * analyzer.dt / 1000.0} ns')
    ave = sum(current[:]) / analyzer.counter
    print(f'Current averaged over full trajectory: {ave}')
    print(f'Norm of average current: {np.linalg.norm(ave)}')

    # Cubic interpolation for computing <J(t)> if there is a constant external field.
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
    figure = plt.figure(figsize=(8, 8))
    plt.subplot(3, 1, 1)
    zeros = np.zeros(shape=analyzer.n_bins)
    plt.plot(t, current_acf, color='red')
    plt.plot(t, zeros, '--', color='black')
    plt.xlabel(r'$t$ (ps)')
    plt.ylabel(r'$\frac{<J(t)J(0)>}{<J(0)J(0)>}$ ((e.nm/ps)$^2$)')

    plt.subplot(3, 1, 2)
    t2 = np.arange(analyzer.dt, (analyzer.counter + 1) * analyzer.dt, step=analyzer.dt)
    # Convert to ns.
    t2 /= 1000.0
    norm_current = np.zeros(shape=(len(current), 1))
    norm_current_t = np.zeros(shape=(len(current), 1))
    zeros2 = np.zeros(shape=norm_current.shape)
    current_t = np.array([0.0, 0.0, 0.0])
    for k in np.arange(0, len(current)):
        current_t += current[k]
        norm_current[k] = np.linalg.norm(current[k])
        ave_current_t = current_t / (k + 1.0)
        norm_current_t[k] = np.linalg.norm(ave_current_t)
    plt.plot(t2, current, color='grey')
    plt.plot(t2, norm_current_t, color='red')
    plt.plot(t2, zeros2, '--', color='black')
    plt.xlabel(r'$t$ (ns)')
    plt.ylabel(r'$J_k, k\in{x,y,z}$')

    if np.linalg.norm(ef_0) > 0.0:
        plt.subplot(3, 1, 3)
        ave = np.zeros(shape=t_new.size)
        a_index = 0
        for a in ave_current:
            ave[a_index] = np.linalg.norm(a)
            a_index += 1
        plt.plot(t_new, ave, color='blue')

    figure.tight_layout(pad=0.5)
    plt.show()
