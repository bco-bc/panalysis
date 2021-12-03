"""Particle system
"""

from particle.particle_spec_catalog import ParticleSpecCatalog
from particle.particle import Particle
import numpy as np


class ParticleSystem():
    """Represents a physical system."""

    def __init__(self):
        self.all = []
        self.free = []
        self.groups = []
        self.box = [0.0, 0.0, 0.0]

    def find(self, pid: str):
        """Finds a particle in this particle system
        :param pid Particle identifier.
        """
        for p in self.all:
            if p.pid == pid:
                return p
        raise ValueError(f'{pid}: No particle with this identifier.')

    def add_free(self, particle: Particle):
        """Adds a free particle to this particle system. The given particle must already exist in this
        particle system
        :param particle Free particle.
        """
        self.find(particle.pid)
        self.free.append(particle)

    def add_particle(self, particle: Particle):
        """Adds another particle to this particle system
        :param particle Particle"""
        self.all.append(particle)

    def add_group(self):
        """Adds another particle group to this particle system
        """
        pass

    def number_of_particles(self) -> int:
        """Returns total number of particles in this particle system
        :return Number
        """
        return len(self.all)


def read_particle_sys(fn: str, catalog: ParticleSpecCatalog) -> ParticleSystem:
    particle_system = ParticleSystem()
    counter = -1
    with open(fn, 'r') as f:
        for line in f:
            items = line.split()
            if counter == -1:
                n_particles = int(items[0])
            if 0 <= counter < n_particles:
                # Read particles
                name = items[0]
                index = int(items[1])
                spec = catalog.find(items[2])
                id = items[3]
                r = np.array([float(items[4]), float(items[5]), float(items[6])])
                v = np.array([float(items[7]), float(items[8]), float(items[9])])
                if spec.protonatable:
                    # Read protonation state.
                    pass
                particle = Particle(id, index, name, spec, r, v)
                particle_system.add_particle(particle)
            # Free particles
            if counter == n_particles:
                n_free = int(items[0])
            if counter == n_particles + 1:
                i_values = np.arange(0, n_free)
                for i in i_values:
                    p = particle_system.find(items[i])
                    particle_system.add_free(p)
            if counter == n_particles + 2:
                n_groups = int(items[0])
            if counter == n_particles + 3 and n_groups != 0:
                # Read groups
                pass
            elif counter == n_particles + 3 and n_groups == 0:
                # Read box
                particle_system.box = np.array([float(items[0]), float(items[1]), float(items[2])])
            counter += 1
    print(f'Particle system holds {len(particle_system.all)} particles.')
    print(f'Particle system holds {len(particle_system.free)} free particles.')
    print(f'Particle system holds {len(particle_system.groups)} particle groups.')
    return particle_system
