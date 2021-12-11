"""Particle system
"""
from typing import List

from particle.particle_spec_catalog import ParticleSpecCatalog
from particle.particle import Particle
from particle.particle_group import ParticleGroup
import numpy as np
import logging


class ParticleSystem:
    """Represents a physical system.
    """

    all: List[Particle]
    free: List[Particle]
    groups: List[ParticleGroup]

    def __init__(self):
        # Lists.
        self.all = []
        self.free = []
        self.groups = []
        self.box = np.zeros(shape=(3, 1))

    def find(self, pid: str):
        """Finds a particle in this particle system
        :param pid Particle identifier.
        """
        for p in self.all:
            if p.pid == pid:
                return p
        raise ValueError(f'{pid}: No particle with this identifier exists in particle system.')

    def add_free(self, particle: Particle):
        """Adds another free particle to this particle system. The given particle must already
        exist in this particle system
        :param particle A free particle.
        """
        self.find(particle.pid)
        self.free.append(particle)

    def add_particle(self, particle: Particle):
        """Adds another particle to this particle system
        :param particle A particle"""
        self.all.append(particle)

    def add_group(self, group: ParticleGroup):
        """Adds another particle group to this particle system
        """
        self.groups.append(group)

    def number_of_particles(self) -> int:
        """Returns total number of particles in this particle system
        :return Number
        """
        return len(self.all)

    def box_volume(self) -> float:
        """Returns box volume
        """
        return self.box[0] * self.box[1] * self.box[2]


def read_particle_sys(fn: str, catalog: ParticleSpecCatalog) -> ParticleSystem:
    """Read particle system from a file
    :param fn File name of file from which particle system is read
    :param catalog Particle catalog
    """

    particle_system = ParticleSystem()
    f = open(fn, 'r')

    # First line
    line = f.readline()
    items = line.split()
    n_particles = int(items[0])
    protonatable = bool(items[1] == '1')

    # Read particles.
    for _ in np.arange(0, n_particles):
        line = f.readline()
        items = line.split()
        name = items[0]
        index = int(items[1])
        spec = catalog.find(items[2])
        pid = items[3]
        r = np.array([float(items[4]), float(items[5]), float(items[6])])
        v = np.array([float(items[7]), float(items[8]), float(items[9])])
        if spec.protonatable:
            # Read protonation state.
            pass
        particle = Particle(pid, index, name, spec, r, v)
        particle_system.add_particle(particle)

    # Read free particles
    line = f.readline()
    items = line.split()
    n_free = int(items[0])
    if n_free > 0:
        line = f.readline()
        items = line.split()
        for k in np.arange(0, n_free):
            p = particle_system.find(items[k])
            particle_system.add_free(p)

    # Read particle groups.
    line = f.readline()
    items = line.split()
    n_groups = int(items[0])
    for _ in np.arange(0, n_groups):
        group = ParticleGroup()

        # Particles in group.
        line = f.readline()
        items = line.split()
        n_particles_in_group = int(items[0])
        line = f.readline()
        items = line.split()
        for k in np.arange(0, n_particles_in_group):
            pid = items[k]
            particle = particle_system.find(pid)
            group.add_particle(particle)
        line = f.readline()
        items = line.split()

        # Bonds.
        n_bonds = int(items[0])
        for _ in np.arange(0, n_bonds):
            line = f.readline()
            items = line.split()
            pid = items[0]
            particle_1 = particle_system.find(pid)
            pid = items[1]
            particle_2 = particle_system.find(pid)
            group.add_bond(p_1=particle_1, p_2=particle_2)

        # Group description compete. Add it to the particle system.
        particle_system.add_group(group)

    # Box
    line = f.readline()
    items = line.split()
    particle_system.box = np.array([float(items[0]), float(items[1]), float(items[2])])

    # Done.
    logging.info(f'Particle system holds {len(particle_system.all)} particles.')
    logging.info(f'Particle system holds {len(particle_system.free)} free particles.')
    logging.info(f'Particle system holds {len(particle_system.groups)} particle groups.')
    logging.info(f'Box dimensions: {particle_system.box}')
    logging.info(f'Protonatable particle system? {protonatable}')

    return particle_system
