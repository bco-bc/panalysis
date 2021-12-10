
from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory
from analysis.current_acf import CurrentACF
import matplotlib.pyplot as plt


if __name__ == '__main__':
    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)

    fn = '/wrk3/tests/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)
    fn = '/wrk3/tests/trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)

    analyzer = CurrentACF(dt=2.0, t_max=500)
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    t, current_acf = analyzer.results()

    plt.plot(t, current_acf, color='red')
    plt.show()
