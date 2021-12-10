"""Plots the total electric current autocorrelation function over given time interval.
"""

from particle.particle_system import read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
from analysis.current_acf import CurrentACF
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate

if __name__ == '__main__':

    if len(sys.argv) < 5:
        print('Plots the total electric current autocorrelation function over given time interval.')
        print('Usage:')
        print('python fnParticleSpecs fnParticleSystem fnTrajectory dt t_max E_0 E_0_direction temperature')
        print('fnParticleSpecs: File name particle specifications')
        print('fnParticleSystem: File name particle system.')
        print('fnTrajectory: File name trajectory file.')
        print('dt: Time interval between trajectory entries')
        print('t_max: Length of time interval.')
        print('E_0: Strength of static homogeneous external electric field E_0')
        print('E_0_direction: Direction of E_0. Either x, y, or z.')
        print('temperature: Temperature of the particle system.')
        raise Exception('Missing argument')

    catalog = read_particle_spec_catalog(sys.argv[1])
    particle_system = read_particle_sys(sys.argv[2], catalog)
    dt = float(sys.argv[4])
    t_max = float(sys.argv[5])
    ef_0_strength = 0.0
    ef_0 = np.zeros(shape=3)
    temperature = 298.15
    if len(sys.argv) > 5:
        ef_0_strength = float(sys.argv[6])
        d = sys.argv[7]
        if d == 'x':
            ef_0[0] = ef_0_strength
        elif d == 'y':
            ef_0[1] = ef_0_strength
        else:
            ef_0[2] = ef_0_strength

    # Compute current ACF.
    trajectory = Trajectory(sys.argv[3])
    analyzer = CurrentACF(dt=dt, t_max=t_max)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    t, current_acf = analyzer.results()
    print(f'Length of time interval: {analyzer.t_max}')
    print(f'Number of states in trajectory: {analyzer.counter}')
    print(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')

    # Cubic interpolation.
    if ef_0_strength > 0.0:
        f_current_acf = interpolate.interp1d(t, current_acf, kind='cubic')

        def current_acf_f(s: float) -> float:
            return f_current_acf(s)

        t_new = np.linspace(0, t[t.size-1], num=2*t.size)
        ave_current = np.zeros(shape=(t_new.size, 3))
        t_index = 0
        val, err = integrate.quadrature(func=current_acf_f, a=0.0, b=t_new[t_new.size-1], maxiter=1000)
        print(f'Value of integration of <J(t)J(0)>: {val}. Error= {err}')
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

    if ef_0_strength > 0.0:
        ave = np.zeros(shape=t_new.size)
        a_index = 0
        for a in ave_current:
            ave[a_index] = np.linalg.norm(a)
            a_index += 1
        plt.plot(t_new, ave, color='blue')

    plt.show()
