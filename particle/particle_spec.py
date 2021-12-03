"""Particle specification
"""


class ParticleSpec():
    """Specifies a particle in terms of mass, charge, etc."""

    def __init__(self,
                 name: str,
                 charge: float,
                 mass: float,
                 radius: float,
                 pKa: float,
                 free: bool,
                 protonatable: bool):
        """Constructor
        :param name Uniquely identifying name
        :param charge Charge value.
        :param mass Mass value.
        :param radius Particle radius.
        :param pKa Acid dissociation constant.
        :param free If true, associated particle is a free particle , i.e. not in any particle
        group
        :param protonatable Particle can bind protons.
        """
        self.name = name
        self.charge = charge
        self.mass = mass
        self.radius = radius
        self.pKa = pKa
        self.free = free
        self.protonatable = protonatable

    def __str__(self):
        pass
