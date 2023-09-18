"""Velocity profile along a particular direction. For instance, along the y-direction in a channel.
"""

import logging

import numpy as np

import util
from analysis.analyzer import Analyzer
from particle.particle_system import ParticleSystem


class VelocityProfile(Analyzer):

    def __init__(self,
                 bin_size: float,
                 box: np.ndarray,
                 direction: str,
                 direction_2: str,
                 pbc: [],
                 exclude=None):
        """Constructor
        :param bin_size Width of each slice
        :param box Simulation box dimensions
        :param direction The direction along which velocity profile should be computed. One of {x,y,z}.
        :param direction_2 The velocity profile is computed perpendicular to this direction. One of
        {x, y, z}
        :param pbc Directions along which PBC should be applied. Subsets of {0,1,2}.
        """
        self.box = box
        self.pbc = pbc
        if direction == 'x':
            self.coordinate = 0
        elif direction == 'y':
            self.coordinate = 1
        else:
            self.coordinate = 2
        if direction_2 == 'x':
            self.coordinate_2 = 0
        elif direction_2 == 'y':
            self.coordinate_2 = 1
        else:
            self.coordinate_2 = 2
        self.n_bins_2 = int(self.box[self.coordinate_2] / bin_size)
        self.bin_size_2 = self.box[self.coordinate_2] / self.n_bins_2
        self.particle_counter = np.zeros(shape=self.n_bins_2)
        self.velocities = np.zeros(shape=self.n_bins_2)
        self.counter = 0
        self.number_of_particles = 0
        self.exclude = exclude

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        if self.counter == 1:
            self.number_of_particles = len(particle_system.all)
        for p in particle_system.all:
            if self.exclude is None or self.exclude != p.spec.name:
                r = util.box.pbc_place_inside(self.box, p.r, pbc=[0, 1, 2])
                index = int(r[self.coordinate_2] / self.bin_size_2)
                if 0 <= index < self.n_bins_2:
                    self.particle_counter[index] += 1
                    v = p.v[self.coordinate]
                    self.velocities[index] += v
                else:
                    logging.warning('Particle position coordinate outside allowed range.')
                    logging.warning(f'Particle ID: {p.pid}')
                    logging.warning(f'Particle position coordinate: {p.r}')
                    logging.warning(f'Index value: {index}')

    def results(self):
        """Returns r, vp(r), the latter is the velocity profile."""
        vp = self.velocities / (self.counter * self.number_of_particles)
        r = np.arange(0, self.box[self.coordinate_2], self.bin_size_2)
        return r, vp
