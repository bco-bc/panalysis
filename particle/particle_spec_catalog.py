"""Particle specifications catalog
"""

from particle.particle_spec import ParticleSpec
import logging


class ParticleSpecCatalog():
    """Holds particle specifications"""

    def __init__(self):
        """Constructor
        """
        self.specs = []

    def add(self, spec: ParticleSpec):
        """Adds another specification to the catalog.
        :param spec Particle specification.
        """
        self.specs.append(spec)

    def find(self, name) -> ParticleSpec:
        """Returns a particle specified. Raises an ValueError if the specification
        cannot be found.
        """
        for spec in self.specs:
            if spec.name == name:
                return spec
        raise ValueError(f"{name}: No such particle specification.")


def read_particle_spec_catalog(fn: str) -> ParticleSpecCatalog:
    """Reads catalog from a file.
    :param fn Filename particle specification catalog.
    """
    catalog = ParticleSpecCatalog()
    counter = 0
    with open(fn, 'r') as f:
        for line in f:
            counter += 1
            if counter > 1:
                items = line.split()
                name = items[0]
                protonatable = items[1] == '1'
                free = items[2] == '1'
                mass = float(items[3])
                charge = float(items[4])
                radius = float(items[5])
                pKa = float(items[6])
                spec = ParticleSpec(name, charge, mass, radius, pKa, free, protonatable)
                catalog.add(spec)

    logging.info(f'Read {counter - 1} particle specifications.')

    return catalog
