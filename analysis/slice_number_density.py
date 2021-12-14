"""Number density in slices
"""

import numpy as np
from analysis.analyzer import Analyzer
from particle.particle_spec import ParticleSpec
from particle.particle_system import ParticleSystem
import util


class SliceNumberDensity(Analyzer):
    """Calculates the number density of particles of a specific particle specification
    in slices parallel to the xy plane, that is in the z-direction"""

    def __init__(self, dz: float, box: np.ndarray, spec: ParticleSpec):
        """Constructor
        :param dz Width of each slice.
        :param box Simulation box dimensions
        :param spec Particle specification for which densities are calculated.
        """
        self.dz = dz
        self.spec = spec
        self.box = box

        self.n = int(self.box[2] / self.dz)
        self.dz = self.box[2] / self.n
        self.histogram = np.zeros(shape=self.n)
        self.counter = 0
        self.nSpec = 0

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        for p in particle_system.all:
            if p.spec.name is self.spec.name:
                r_in = util.box.pbc_place_inside(particle_system.box, p.r)
                z = r_in[2]
                index = int(z / self.dz)
                self.histogram[index] += 1.0

    def results(self):
        """Returns number density """
        volume_slice = self.box[0] * selolef.box[1] * self.dz
        number_densities = self.histogram / (float(self.counter) * volume_slice)
        z = np.zeros(shape = self.n)
        k_values = np.arange(0, self.n)
        for k in k_values:
            z[k] = k * self.dz
        return z, number_densities

