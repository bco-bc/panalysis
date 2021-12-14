import numpy as np
from particle.particle_spec import ParticleSpec


class Particle:
    """Holds a particle"""

    def __init__(self,
                 pid: str,
                 index: int,
                 name: str,
                 spec: ParticleSpec,
                 r: np.ndarray,
                 v: np.ndarray):
        """Constructor
        :param pid Particle unique identifier
        :param index Particle index >= 0.
        :param name Particle's name.
        :param spec Particle specification.
        :param r Position.
        :param v Velocity.
        """
        self.pid = pid
        self.index = index
        self.name = name
        self.spec = spec
        self.r = r
        self.v = v

    def charge(self) -> float:
        """Returns particle's charge
        """
        return self.spec.charge

    def mass(self) -> float:
        """Returns particle's mass
        """
        return self.spec.mass
