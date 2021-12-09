"""Plot a particle's displacement over a given time interval.
"""

from particle.particle_system import read_particle_sys
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
from analysis.msd import MSD
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


if __name__ == '__main__':

    if len(sys.argv) < 6:
        print('Usage:')
        print('python fnParticleSpecs fnParticleSystem fnTrajectory spec dt t_max')
        print('fnParticleSpecs: File name particle specifications')
        print('fnParticleSystem: File name particle system.')
        print('fnTrajectory: File name trajectory file.')
        print('spec: Particle specification name.')
        print('dt: Time interval between trajectory entries')
        print('t_max: Length of time interval for displacement calculation.')
        raise Exception('Missing argument')

    catalog = read_particle_spec_catalog(sys.argv[1])
    particle_system = read_particle_sys(sys.argv[2], catalog)
    spec = catalog.find(sys.argv[4])
    dt = float(sys.argv[5])
    t_max = float(sys.argv[6])

    # Compute displacement.
    trajectory = Trajectory(sys.argv[3])
    analyzer = MSD(dt, t_max, spec)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    t, msd = analyzer.results()
    print(f'Number of particles with particle specification \'{analyzer.spec.name}\': {analyzer.n_specs}')
    print(f'Length of time interval: {analyzer.t_max}')
    print(f'Number of states in trajectory: {analyzer.counter}')
    print(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')

    x = t.reshape((-1, 1))
    model = LinearRegression()
    model.fit(x, msd)
    r_sq = model.score(x, msd)
    print("Coefficients of determination:", r_sq)
    a = model.coef_[0]
    b = model.intercept_
    print("Intercept: ", b)
    print("Slope: ", a)
    prediction = np.zeros(shape=len(t))
    index = 0
    for v in t:
        prediction[index] = a * v + b
        index += 1

    D = prediction[t.size - 1] / (6.0 * t[t.size - 1])
    print(f'Diffusion (self) constant: {D} nm^2/ps')

    # Plot displacement.
    plt.plot(t, msd, color='red')
    plt.plot(t, prediction, color='blue')
    plt.xlabel(r'$t$ (ps)')
    plt.ylabel(r'$<|r(t+s)r(s)|^2>$ (nm$^2$)')
    plt.show()
