"""Plots the total electric current (flux) autocorrelation function over given time interval.
"""

from particle.particle_system import read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import ParticleSystem
from simulation.trajectory import Trajectory
from analysis.current_acf import CurrentACF
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate, constants
import logging
import util


def usage():
    print()
    print('Usage: python plot_current_acf.py [OPTION [value]]...')
    print('Description: Plots total electric current (flux) autocorrelation function.')
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
        d = conf['direction']
        if d == 'x':
            ef_0[0] = e_0_strength
        elif d == 'y':
            ef_0[1] = e_0_strength
        else:
            ef_0[2] = e_0_strength
    except (KeyError, FileExistsError, FileNotFoundError):
        usage()
        sys.exit("Missing input argument(s).")

    # Perform the analysis.
    trajectory = Trajectory(fn_trajectory)
    analyzer = CurrentACF(dt=dt, t_max=t_max)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    # Get the results.
    t_current_current, t_current_acf, cos_a_probability_density = analyzer.results()
    logging.info(f'Length of time interval: {analyzer.t_max}')
    logging.info(f'Number of states/entries in trajectory: {analyzer.counter}')
    logging.info(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')
    logging.info(f'Length of trajectory: {analyzer.counter * analyzer.dt / 1000.0} ns')

    t, current = t_current_current
    ave = sum(current[:]) / analyzer.counter
    print(f'Current averaged over full trajectory: {ave}')
    print(f'Norm of average current: {np.linalg.norm(ave)}')

    # Cubic interpolation for computing <J(t)> and the conductivity.
    t_acf, current_acf = t_current_acf
    f_current_acf = interpolate.interp1d(t_acf, current_acf, kind='cubic')
    n = t_acf.size * 2
    t_new = np.linspace(t_acf[0], t_acf[t_acf.size-1], num=n, endpoint=True)
    acf_new = f_current_acf(t_new)

    def current_acf_f(s) -> float:
        """ Function for integration of ACF
        :param s:
        :return: Value of ACF at time s.
        """
        return f_current_acf(s)

    ave_current = np.zeros(shape=(t_new.size, 3))
    t_index = 0
    a = 0
    b = t_acf[t_acf.size-1]
    logging.info(f'Integrating ACF from {a} to {b}.')
    val, err = integrate.quadrature(func=current_acf_f, a=a, b=b, maxiter=1000)
    logging.info(f'Value of integration of <J(t)J(0)>: {val}. Error= {err}')
    k_B = constants.k * constants.N_A / 1000.0  # To kJ/(mol K)
    conductivity = val / (3.0 * k_B * temperature * particle_system.box_volume())
    print(f'Conductivity: {conductivity}')

    # Plot results.
    figure = plt.figure(figsize=(8, 8))

    # Current ACF
    plt.subplot(4, 1, 1)
    zeros = np.zeros(shape=analyzer.n_bins)
    plt.plot(t_acf, current_acf, 'o', color='red')
    plt.plot(t_acf, zeros, '--', color='black')
    plt.plot(t_new, acf_new, '--', color='blue')
    plt.xlabel(r'$t$ (ps)')
    plt.ylabel(r'$\frac{<J(t)J(0)>}{<J(0)J(0)>}$')

    # Current
    plt.subplot(4, 1, 2)
    t /= 1000.0  # Conversion to ns.
    norm_current = np.zeros(shape=(len(current), 1))
    norm_current_t = np.zeros(shape=(len(current), 1))
    zeros2 = np.zeros(shape=norm_current.shape)
    current_t_sum = np.array([0.0, 0.0, 0.0])
    for k in np.arange(0, len(current)):
        current_t_sum += current[k]
        norm_current[k] = np.linalg.norm(current[k])
        ave_current_t = current_t_sum / (k + 1.0)
        norm_current_t[k] = np.linalg.norm(ave_current_t)
    plt.plot(t, current, color='grey')
    plt.plot(t, norm_current_t, color='red')
    plt.plot(t, zeros2, '--', color='black')
    plt.xlabel(r'$t$ (ns)')
    plt.ylabel(r'$J_k, k\in{x,y,z}$')

    if np.linalg.norm(ef_0) > 0.0:
        # Plot <J(t)> as function of t
        # plt.subplot(4, 1, 3)
        pass

    # Probability density function of cos(angle).
    plt.subplot(4, 1, 4)
    cos_a, probability_density = cos_a_probability_density
    plt.plot(cos_a, probability_density, color='blue')
    plt.xlabel(r'$\cos(\alpha)$')
    plt.ylabel(r'$p(\cos(\alpha))$')
    plt.ylim(0.0, 1.0)

    figure.tight_layout(pad=1.0)
    plt.show()
