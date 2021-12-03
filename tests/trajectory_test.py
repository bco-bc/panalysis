from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory


if __name__ == '__main__':
    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)
    fn = '/home/andre/simulations/electrolyte/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)

    fn = '/home/andre/simulations/electrolyte/trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)
    counter = 0
    while trajectory.next(particle_system):
        counter += 1
        print(f'counter = {counter}')
    trajectory.close()

