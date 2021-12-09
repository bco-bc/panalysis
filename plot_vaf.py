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


if __name__ == '__main__':

    if len(sys.argv) < 7:
        print('Usage:')
        print('python fnParticleSpecs fnParticleSystem fnTrajectory spec dt t_max temperature')
        print('fnParticleSpecs: File name particle specifications')
        print('fnParticleSystem: File name particle system.')
        print('fnTrajectory: File name trajectory file.')
        print('spec: Particle specification name.')
        print('dt: Time interval between trajectory entries')
        print('t_max: Length of time interval for displacement calculation.')
        print('temperature: Temperature of the particle system. Default is 298.18')
        raise Exception('Missing argument')

    catalog = read_particle_spec_catalog(sys.argv[1])
    particle_system = read_particle_sys(sys.argv[2], catalog)
    spec = catalog.find(sys.argv[4])
    dt = float(sys.argv[5])
    t_max = float(sys.argv[6])
    temperature = 298.15
    if len(sys.argv) == 8:
        temperature = float(sys.argv[7])

    # Compute VAF.
    trajectory = Trajectory(sys.argv[3])
    bin_size = dt
    analyzer = VAF(dt, t_max, spec)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    t, vaf = analyzer.results()
    print(f'Number of particles with particle specification \'{analyzer.spec.name}\': {analyzer.n_specs}')
    print(f'Length of time interval: {analyzer.t_max}')
    print(f'Number of states in trajectory: {analyzer.counter}')
    print(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')

    # Create cubic interpolation.
    f_vaf = interpolate.interp1d(t, vaf, kind='cubic')
    n = t.size * 2
    t_new = np.linspace(0, t[t.size - 1], num=2*t.size, endpoint=True)
    vaf_new = f_vaf(t_new)

    def vaf_f(time) -> float:
        return f_vaf(time)

    a = 0.0
    b = t[t.size-1]
    val, err = integrate.quadrature(func=vaf_f, a=a, b=b, maxiter=1000)
    k = constants.k * constants.N_A / 1.0e+03  # In kJ/(mol K)
    D = k * temperature / spec.mass * val
    print(f'Diffusion (self) constant: {D} nm^2/ps')

    # Plot VAF
    zeros = np.zeros(shape=analyzer.n_bins)
    plt.plot(t, vaf, 'o', color='red')
    plt.plot(t_new, vaf_new, '--', color='blue')
    plt.xlabel(r'$t$ (ps)')
    plt.plot(t, zeros, color='black')
    plt.ylabel(r'$<v(t+s)v(t)$')
    plt.legend(['data', 'cubic'])
    plt.show()
