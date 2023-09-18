import numpy as np

from particle.particle_spec import ParticleSpec
from particle.particle_system import ParticleSystem


class Analyzer:

    def perform(self, particle_system: ParticleSystem):
        """Perform an analysis of the given particle system
        :param particle_system Particle system
        """
        pass

    def results(self):
        """Return the results
        """
        pass

    @staticmethod
    def number_of_specs(particle_system: ParticleSystem, spec: ParticleSpec) -> int:
        """Returns number of particle of given specification.
        """
        n_specs = 0
        for p in particle_system.all:
            if p.spec.name == spec.name:
                n_specs += 1
        return n_specs

    @staticmethod
    def positions(particle_system: ParticleSystem) -> np.ndarray:
        """Returns current positions
        """
        size = len(particle_system.all)
        r = np.ndarray(shape=(size, 3))
        for i in np.arange(0, size):
            r[i] = particle_system.all[i].r
        return r

    @staticmethod
    def velocities(particle_system: ParticleSystem) -> np.ndarray:
        """Returns current velocities.
        """
        size = len(particle_system.all)
        v = np.ndarray(shape=(size, 3))
        for i in np.arange(0, size):
            v[i] = particle_system.all[i].v
        return v
