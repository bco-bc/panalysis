from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys


if __name__ == '__main__':
    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)
    fn = '/wrk3/tests/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)
    print(f'Number of particles: {particle_system.number_of_particles()}')
