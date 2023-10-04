"""Distance between specific particles
"""

from analysis.analyzer import Analyzer
from particle.particle_system import ParticleSystem
import util
import numpy as np


class Distance(Analyzer):
    """Monitors the distance between two specific particles or between particles of given specification. For the
    latter, only pairs within the -same- particle group are considered.
    """

    def __init__(self, dt: float = 0.002,
                 id_1: str = '',
                 id_2: str = '',
                 pbc: [] = None,
                 spec_1: str = '',
                 spec_2: str = ''):
        """Constructor.
        :param dt: Time difference between states
        :param id_1: Identifier particle #1. If provided, then id_2 must be provided as well, and spec_1 and spec_2
        must not.
        :param id_2: Identifier particle #2. If provided, then id_1 must be provided as well, and spec_1 and spec_2
        must not.
        :param pbc: Directions along which PBC should be applied. Subsets of {0,1,2}
        :param spec_1: Specification #1. If provided, then spec_2 must be provided as well, and id_1 and id_2
        must not.
        :param spec_2: Specification #2 If provided, then spec_1 must be provided as well, and id_1 and id_2
        must not.
        """
        self.dt = dt
        self.id_1 = id_1
        self.id_2 = id_2
        self.index_1 = 0
        self.index_2 = 0
        if pbc is None:
            self.pbc = [0, 0, 0]
        else:
            self.pbc = pbc
        self.spec_1 = spec_1
        self.spec_2 = spec_2

        # Two specific particles? Or just all particle pairs of particular specification within particle
        # groups?
        self.specs = len(self.spec_1) and len(self.spec_2)
        print(f'{self.specs}: Specification-based distance?')
        if self.specs:
            print(f'{self.spec_1}: Particle specification #1')
            print(f'{self.spec_2}: Particle specification #2')
        else:
            print(f'{self.id_1}: Particle identifier #1')
            print(f'{self.id_2}: Particle identifier #2')

        # Set up.
        self.counter = 0
        self.distances: list = []
        self.x: list = []
        self.y: list = []
        self.z: list = []

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1

        if self.counter == 1:
            if not self.specs:
                pi = particle_system.find(self.id_1)
                pj = particle_system.find(self.id_2)
                self.index_1 = pi.index
                self.index_2 = pj.index

        if not self.specs:
            particles = particle_system.all
            pi = particles[self.index_1]
            pj = particles[self.index_2]
            r1 = pi.r
            r2 = pj.r
            r_ij = util.box.pbc_distance(particle_system.box, r1, r2, self.pbc)
            distance = np.linalg.norm(r_ij)
            self.distances.append(distance)
            self.x.append(r_ij[0])
            self.y.append(r_ij[1])
            self.z.append(r_ij[2])
        else:
            for group in particle_system.groups:
                for i in np.arange(0, len(group.particles) - 1):
                    pi = group.particles[i]
                    if pi.spec.name == self.spec_1:
                        ri = pi.r
                        for j in np.arange(i + 1, len(group.particles)):
                            pj = group.particles[j]
                            if pj.spec.name == self.spec_2:
                                rj = pj.r
                                r_ij = util.box.pbc_distance(particle_system.box, ri, rj, self.pbc)
                                distance = np.linalg.norm(r_ij)
                                self.distances.append(distance)
                                self.x.append(r_ij[0])
                                self.y.append(r_ij[1])
                                self.z.append(r_ij[2])

    def results(self) -> (np.ndarray, np.ndarray):
        """Returns t, distance[t], x(t), y(t), z(t) when dealing specific particles. Otherwise,
        distances, x, y, z, are returned
        """
        if not self.specs:
            t = np.zeros(shape=len(self.distances))
            for k in np.arange(0, len(self.distances)):
                t[k] = k * self.dt
            return t, self.distances, self.x, self.y, self.z
        else:
            return self.distances, self.x, self.y, self.z
