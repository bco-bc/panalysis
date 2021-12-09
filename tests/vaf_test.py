from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory
from analysis.vaf import VAF
import matplotlib.pyplot as plt

if __name__ == '__main__':
    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)

    fn = '/home/juffer/simulations/electrolyte/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)
    fn = '/home/juffer/simulations/electrolyte/trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)

    analyzer = VAF(dt=20.0, t_max=500.0, spec=catalog.find('Na+'))
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    t, vaf = analyzer.results()

    plt.plot(t, vaf, color='red')
    plt.show()
