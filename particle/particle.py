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
        self.pid = pid
        self.index = index
        self.name = name
        self.spec = spec
        self.r = r
        self.v = v
