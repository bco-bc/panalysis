"""Particle group
"""
from typing import List,Set
import numpy as np

import util.box
from particle.particle import Particle


class ParticleGroup:
    """Group of particles forming a logical unit.
    """
    particles: List[Particle]
    bonds: List[Set[Particle]]

    def __init__(self):
        """Constructor"""
        self.particles = []
        self.bonds = []

    def add_particle(self, particle: Particle):
        """Adds another particle to this group
        :param particle A particle
        """
        self.particles.append(particle)

    def add_bond(self, p_1: Particle, p_2: Particle):
        """Adds a bond between two particles to this group
        :param p_1 First particle
        :param p_2 Second particle
        """
        self.bonds.append({p_1, p_2})

    def find(self, pid: str) -> Particle:
        """Finds particle with given particle identifier."""
        for p in self.particles:
            if p.pid == pid:
                return p
        raise ValueError(f'{pid}: No particle with this identifier exists in particle group.')

    def dipole_moment(self, box: np.ndarray, pbc:[]) -> np.ndarray:
        """Returns the dipole moment of this group relative to the origin.
        :param box Simulation box.
        :param pbc  :param pbc Directions specification. Subsets of {0,1,2}.
        """
        m_total = np.zeros(shape=3)
        origin = np.zeros(shape=3)
        for p in self.particles:
            r = util.box.pbc_distance(box, p.r, origin, pbc)
            m_total += p.charge() * r
        return m_total

    def position(self):
        """Returns the center of mass.
        """
        total_mass = 0.0
        cm = np.zeros(shape=3)
        for p in self.particles:
            cm +=  p.mass() * p.r
            total_mass += p.mass()
        return cm / total_mass
