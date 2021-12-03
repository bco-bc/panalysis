from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory
from analysis.slice_number_density import SliceNumberDensity


if __name__ == '__main__':
    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)
    fn = '/home/andre/simulations/electrolyte/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)

    fn = '/home/andre/simulations/electrolyte/trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)

    analyzer = SliceNumberDensity(0.1, particle_system.box, catalog.find('Na+'))
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    print(f'Result: {analyzer.results()}')
    
