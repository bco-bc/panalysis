from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import logging

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG)

    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)

    # fn = '/wrk3/tests/electrolyte.ps'
    fn = '/localdisk/resources/small.ps'
    particle_system = read_particle_sys(fn, catalog)
