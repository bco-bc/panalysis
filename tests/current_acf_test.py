import numpy as np

from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory
from analysis.current_acf import CurrentACF
import matplotlib.pyplot as plt

if __name__ == '__main__':
    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)

    fn = '/home/juffer/simulations/electrolyte/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)
    fn = '/home/juffer/simulations/electrolyte//trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)

    analyzer = CurrentACF(dt=2.0, t_max=500)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    t, current, current_acf = analyzer.results()
    ave = sum(current[:]) / analyzer.counter
    print(f'Average current: {ave}')
    print(f'Norm of average current: {np.linalg.norm(ave)}')

    figure = plt.figure(figsize=(8, 8))

    plt.subplot(2, 1, 1)
    plt.plot(t, current_acf, color='red')

    plt.subplot(2, 1, 2)
    t2 = np.arange(analyzer.dt, (analyzer.counter + 1) * analyzer.dt, step=analyzer.dt)
    norm_current = np.zeros(shape=(len(current), 1))
    for k in np.arange(0, len(current)):
        norm_current[k] = np.linalg.norm(current[k])
    # plt.plot(t2, current, color='red')
    plt.plot(t2, norm_current, color='blue')

    plt.show()
