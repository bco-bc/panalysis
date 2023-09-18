"""Plots total dipole moment autocorrelation function over given time interval.
"""

import math
import logging
import sys
import util
from simulation.trajectory import Trajectory
from particle.particle_system import read_particle_sys, ParticleSystem
from particle.particle_spec_catalog import read_particle_spec_catalog
from analysis.dipole_moment import DipoleMoment
from analysis import analyze
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq

def usage():
    print()
    print('Usage: python plot_dipole_moment.py [OPTION [value]]...')
    print('Description: Plots total dipole moment autocorrelation function.')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('-dt VALUE: VALUE is the time difference between states in trajectory.')
    print('-t-max VALUE: VALUE is the length of time interval over which dipole '
          'moment auto correlation function is calculated.')
    print()
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print('-T VALUE: VALUE is the temperature. Default is 298.15.')
    print('-n-skip VALUE: Skip the first VALUE state(s) from trajectory. Default is 0.')
    print('-bin-size VALUE: Value is the bin size for the theta probability density function. '
          'Default is 0.1')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)
    conf = util.parse_argv.parse(sys.argv)
    logging.info(f'Received program arguments {conf}')
    fn_trajectory = 'trajectory.dat'
    t_max = 0.0
    spec = ''
    dt = 0.0
    particle_system = ParticleSystem()
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        dt = conf['dt']
        t_max = conf['t-max']
        fn_trajectory = conf['fn-trajectory']
        pbc = conf['pbc']
        n_skip = int(conf['n-skip'])
        bin_size = float(conf['bin-size'])
    except (KeyError, FileExistsError, FileNotFoundError):
        usage()
        sys.exit("Missing input argument(s) and/or non-existing files.")

    # Perform analysis.
    trajectory = Trajectory(fn_trajectory)
    analyzer = DipoleMoment(dt=dt, t_max=t_max, pbc=pbc, bin_size=bin_size)
    analyze.perform(analyzer, particle_system, trajectory, n_skip)
    trajectory.close()

    # Get results
    t_m, ave_m, thetas, cos_thetas, theta_density, cos_theta_density, t_m_acf, n_states = analyzer.results()
    print(f'Number of states employed for analysis: {n_states}')

    # Average dipole moment.
    print(f'Dipole moment averaged of trajectory: {ave_m}')
    print(f'Norm average dipole moment: {np.linalg.norm(ave_m)}')
    print(f'Polar (theta) angle average dipole moment: {util.cvector.polar_angle(ave_m)}')
    print(f'cos of polar angle average dipole moment: {util.cvector.cos_polar_angle(ave_m)}')

    # Plot results.
    figure = plt.figure(figsize=(8, 8))

    # Dipole moment ACF.
    plt.subplot(4, 2, 1)
    t_acf, m_acf = t_m_acf
    zeros = np.zeros(shape=analyzer.n_bins)
    plt.plot(t_acf, m_acf, color='red')
    plt.plot(t_acf, zeros, '--', color='black')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\frac{<M(t)M(0)>}{<M(0)M(0)>}$')

    # Dipole moment.
    plt.subplot(4, 2, 2)
    t, dipole_moment = t_m
    norm_m = np.zeros(shape=len(dipole_moment))
    norm_m_t = np.zeros(shape=len(dipole_moment))
    m_t = np.array([0.0, 0.0, 0.0])
    m_x = np.zeros(shape=len(dipole_moment))
    m_y = np.zeros(shape=len(dipole_moment))
    m_z = np.zeros(shape=len(dipole_moment))
    for i in np.arange(0, len(dipole_moment)):
        m_t += dipole_moment[i]
        norm_m[i] = np.linalg.norm(dipole_moment[i])
        ave_m_t = m_t / (1.0 + i)
        norm_m_t[i] = np.linalg.norm(ave_m_t)
        m_x[i] = dipole_moment[i][0]
        m_y[i] = dipole_moment[i][1]
        m_z[i] = dipole_moment[i][2]
    zeros = np.zeros(shape=t.shape)
    # plt.plot(t, dipole_moment, color='grey')
    plt.plot(t, m_x, color='orange', label=r'M$_x$')
    plt.plot(t, m_y, color='blue', label=r'M$_y$')
    plt.plot(t, m_z, color='green', label=r'M$_z$')
    plt.plot(t, norm_m_t, color='red')
    plt.plot(t, zeros, '--', color='black')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$M(t)_k, k \in{x,y,z}$')
    plt.legend()

    # Probability density function of theta.
    plt.subplot(4, 2, 3)
    theta, p_theta = theta_density
    halve_sin_theta = np.zeros(len(theta))
    for k in np.arange(0, len(theta)):
        halve_sin_theta[k] = 0.5 * math.sin(theta[k])
    plt.plot(theta, p_theta, color='blue')
    plt.plot(theta, halve_sin_theta, color='red')
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$p(\theta)$')

    # Probability density function of cos(theta).
    plt.subplot(4, 2, 4)
    cos_theta, p_cos_theta = cos_theta_density
    halve = np.full(shape=len(cos_theta), fill_value=0.5)
    plt.plot(cos_theta, p_cos_theta, color='blue')
    plt.plot(cos_theta, halve, color='red')
    plt.ylabel(r'$p(\cos(\theta))$')
    plt.xlabel(r'$\cos(\theta)$')

    # Thetas versus time
    plt.subplot(4, 2, 5)
    t = np.zeros(len(thetas))
    for k in np.arange(0, len(thetas)):
        t[k] = k * dt
    plt.plot(t, thetas, color='red')
    plt.xlabel(r't')
    plt.ylabel(r'$\theta$')
    average = sum(thetas[:]) / len(thetas)
    print(f'Average value of theta from theta data set: {average}')
    print(f'Value of pi/2: {0.5 * math.pi}')

    # Fourier transform of ACF
    plt.subplot(4, 2, 6)
    N = len(m_acf)
    T = t_acf[1] - t_acf[0]
    x = t_acf
    y = m_acf
    yf = fft(y)
    xf = fftfreq(N, T)[:N//2]
    plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
    #plt.xlim(0,2)
    plt.xlabel('f')
    plt.ylabel('F')



    # print(f'Frequencies: {xf}')
    # print(f'Transform: {yf}')

    # Simple integration of probability density functions.
    integral = 0.0
    for k in np.arange(len(p_theta) - 1):
        integral += 0.5 * (p_theta[k] + p_theta[k+1]) * bin_size
    average = 0.0
    for k in np.arange(len(p_theta)):
        average += theta[k] * bin_size * p_theta[k]
    print(f'Integral of probability density function p(theta) (must be 1): {integral}')
    print(f'Average value of theta from probability density function p(theta): {average}')

    integral = 0.0
    average = 0.0
    for k in np.arange(len(p_cos_theta) - 1):
        cos_bin_size = math.fabs(cos_theta[k] - cos_theta[k+1])
        integral += 0.5 * (p_cos_theta[k] + p_cos_theta[k+1]) * cos_bin_size
        average += cos_theta[k] * p_cos_theta[k] * cos_bin_size
    print(f'Integral of probability density function p(cos(theta)) (must be 1): {integral}')
    print(f'Average value of cos(theta) from probability density function p(cos(theta)): {average}')

    # print()
    # print(f'theta: {thetas}')
    # print(f'p_theta: {p_theta}')
    # print()
    # print(f'cos_theta: {cos_thetas}')
    # print(f'p_cos_theta: {p_cos_theta}')

    figure.tight_layout(pad=1.0)
    plt.show()
