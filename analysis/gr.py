"""Radial distribution function
"""

import math
import numpy as np
from analysis.analyzer import Analyzer
from particle.particle_system import ParticleSystem
from particle.particle_spec import ParticleSpec
import util

class Gr(Analyzer):
    """Calculates the radial distribution g(r) between two particles of given
    specification."""

    def __init__(self,
                 bin_size: float,
                 box: np.ndarray,
                 r_max: float,
                 spec_1: ParticleSpec,
                 spec_2: ParticleSpec):
        """Constructor. All argument are required.
        :param bin_size Bin size ofor g(r)
        :param box Simulation box dimensions.
        :param r_max Maximum r for g(r). Should be smaller than halve of any box dimensions.
        :param spec_1 First of two particle specifications.
        :param spec_2 Second of two particle specifications.
        """
        self.bin_size = bin_size
        self.box = box
        self.r_max = r_max
        self.spec_1 = spec_1
        self.spec_2 = spec_2

        self.n_bins = int(self.r_max / self.bin_size)
        self.bin_size = self.r_max / self.n_bins
        self.histogram = np.zeros(shape=self.n_bins)
        self.counter = 0
        self.n_spec_1 = 0
        self.n_spec_2 = 0

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        if self.counter == 1:
            for p in particle_system.all:
                if p.spec.name is self.spec_1.name:
                    self.n_spec_1 += 1
                if p.spec.name is self.spec_2.name:
                    self.n_spec_2 += 1
            if self.n_spec_1 == 0 or self.n_spec_2 == 0:
                raise ValueError(f'({self.spec_1.name}, {self.spec_2.name}): no such particle(s).')
        for pi in particle_system.all:
            if pi.spec.name == self.spec_1.name:
                r_i = pi.r
                for pj in particle_system.all:
                    if pi.pid != pj.pid:
                        if pj.spec.name == self.spec_2.name:
                            r_j = pj.r
                            r_ij = util.box.pbc_distance(self.box, r_i, r_j)
                            distance = np.linalg.norm(r_ij)
                            index = int(distance / self.bin_size)
                            if index < self.histogram.size:
                                self.histogram[index] += 1

    def results(self) -> (np.ndarray, np.ndarray):
        """Returns r and g(r)
        """
        gr = np.zeros(shape = self.n_bins)
        r = np.zeros(shape = self.n_bins)
        factor = 4.0 * math.pi / 3.0
        rho_2 = self.n_spec_2 / (self.box[0] * self.box[1] * self.box[2])
        i_values = np.arange(0, self.n_bins)
        for i in i_values:
            r_i = i * self.bin_size
            r[i] = r_i
            r_ii = (i + 1) * self.bin_size
            d_v = factor * (r_ii * r_ii * r_ii - r_i * r_i * r_i)
            n_2 = rho_2 * d_v
            if self.counter > 0 and self.n_spec_2 > 0:
                gr[i] = float(self.histogram[i]) / float(n_2 * self.n_spec_1 * self.counter)
            else:
                gr[i] = 0.0
        return r, gr
