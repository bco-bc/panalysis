"""Probability density function (pdf) of particles in slices.
"""

import numpy as np
from analysis.analyzer import Analyzer
from particle.particle_spec import ParticleSpec
from particle.particle_system import ParticleSystem
import logging
import util


class SliceNumberDensity(Analyzer):
    """Calculates the probability of observing a particle of a specific particle specification
    in slices parallel to the xy/yz/zx plane, that is in the z/y/z-direction"""

    def __init__(self,
                 bin_size: float,
                 box: np.ndarray,
                 r_max: float,
                 spec: ParticleSpec,
                 direction: str,
                 pbc: []):
        """Constructor
        :param bin_size Width of each slice
        :param box Simulation box dimensions
        :param r_max Maximum r for probability density function (pdf). Should not be larger than the box
         dimension along which the pdf is computed
        :param spec Particle specification for which densities are calculated
        :param direction The direction along which pdf is computed. One of {x,y,z}.
        :param pbc Directions along which PBC should be applied. Subsets of {0,1,2}.
        """
        self.bin_size = bin_size
        self.box = box
        self.spec = spec
        self.r_max = r_max
        if direction == 'x':
            self.coordinate = 0
        elif direction == 'y':
            self.coordinate = 1
        else:
            self.coordinate = 2
        self.pbc = pbc

        self.n_bins = int(self.r_max / self.bin_size)
        self.bin_size = self.r_max / self.n_bins
        self.histogram = np.zeros(shape=self.n_bins)
        self.counter = 0
        self.nSpec = 0

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        for p in particle_system.all:
            if p.spec.name is self.spec.name:
                r = util.box.pbc_place_inside(self.box, p.r, pbc=[0, 1, 2])
                index = int(r[self.coordinate] / self.bin_size)
                if 0 <= index < self.n_bins:
                    self.histogram[index] += 1.0
                else:
                    logging.warning('Particle position outside allowed range.')
                    logging.warning(f'Particle ID: {p.pid}')
                    logging.warning(f'Particle position: {p.r}')
                    logging.warning(f'Index value: {index}')

    def results(self):
        """Returns probability density function (pdf)
        """
        pdf = self.histogram / self.counter
        total = sum(pdf[:])
        pdf_normalized = pdf / (total * self.bin_size)
        r = np.zeros(shape=self.n_bins)
        for k in np.arange(0, self.n_bins):
            r[k] = k * self.bin_size
        return r, pdf_normalized
