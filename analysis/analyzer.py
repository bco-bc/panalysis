import numpy as np

from particle.particle_spec import ParticleSpec
from particle.particle_system import ParticleSystem


class Analyzer:

    def perform(self, particle_system: ParticleSystem):
        """Perform an analysis of the given particle system
        :param particle_system Particle system
        """
        pass

    @staticmethod
    def number_of_specs(particle_system: ParticleSystem, spec: ParticleSpec):
        """Returns number of particle of given specification.
        """
        n_specs = 0
        for p in particle_system.all:
            if p.spec.name == spec.name:
                n_specs += 1
        return n_specs

    @staticmethod
    def positions(particle_system: ParticleSystem):
        """Returns current positions
        """
        r = np.ndarray(shape=(len(particle_system.all), 3))
        p_index = 0
        for p in particle_system.all:
            r[p_index] = p.r
            p_index += 1
        return r

    @staticmethod
    def velocities(particle_system: ParticleSystem) -> np.ndarray:
        """Returns current velocities.
        """
        v = np.ndarray(shape=(len(particle_system.all), 3))
        p_index = 0
        for p in particle_system.all:
            v[p_index] = p.v
            p_index += 1
        return v
