"""Plot a particle's displacement over a given time interval.
"""

from particle.particle_system import read_particle_sys, ParticleSystem
from particle.particle_spec_catalog import read_particle_spec_catalog
from simulation.trajectory import Trajectory
from analysis.msd import MSD
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import util
import logging


def usage():
    print()
    print('Usage: python plot_msd.py [OPTION [value]]...')
    print('Description: Computes and plots the mean square displacement msd(t) over '
          'given time interval.')
    print()
    print('Required arguments:')
    print('-ps FN: FN is the file name of the particle system')
    print('-dt VALUE: VALUE is the time difference (in ps) between states in trajectory.')
    print('-t-max VALUE: VALUE is the length of time interval (in ps) over which VAF is calculated.')
    print('-spec NAME: NAME is the particle specification names.')
    print()
    print('Optional arguments:')
    print('-s FN: FN is the filename particle specifications. Default is \'particle-specs.dat\'.')
    print('-tr FN: FN is the trajectory file name. Default is \'trajectory.dat\'')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    logging.info(f'Received program arguments {conf}')
    fn_trajectory = 'trajectory.dat'
    temperature = 298.15
    t_max = 0.0
    dt = 0.0
    particle_system = ParticleSystem()
    spec = None
    try:
        catalog = read_particle_spec_catalog(conf['fn-particle-specs'])
        fn_ps = conf['fn-particle-system']
        particle_system = read_particle_sys(fn_ps, catalog)
        dt = conf['dt']
        t_max = conf['t-max']
        temperature = conf['temperature']
        fn_trajectory = conf['fn-trajectory']
        spec = catalog.find(conf['spec'])
    except (KeyError, FileExistsError, FileNotFoundError):
        usage()
        sys.exit("Missing input argument(s).")

    # Compute displacement.
    trajectory = Trajectory(fn_trajectory)
    analyzer = MSD(dt, t_max, spec)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()
    t, msd = analyzer.results()
    print(f'Results: Times: {t}, MSD: {msd}')
    logging.info(f'Number of particles with particle specification \'{analyzer.spec.name}\': {analyzer.n_specs}')
    logging.info(f'Length of time interval: {analyzer.t_max}')
    logging.info(f'Number of states in trajectory: {analyzer.counter}')
    logging.info(f'Length of trajectory: {analyzer.counter * analyzer.dt} ps')

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
    print(f'(Self) Diffusion constant: {D} nm^2/ps')
    D *= 1.0e-02
    print(f'(Self) Diffusion constant: {D} cm^2/s')

    # Plot displacement.
    plt.plot(t, msd, color='red')
    plt.plot(t, prediction, color='blue')
    plt.xlabel(r'$t$ (ps)')
    plt.ylabel(r'$<|r(t+s)r(s)|^2>$ (nm$^2$)')
    plt.show()
