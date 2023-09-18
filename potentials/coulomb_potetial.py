"""Regular Coulomb interaction between two point charges.
U(r) = Q1 * Q2 / (4.0 * pi * eps0 * r), where Q1 and Q2 are charge values of the two point charges, r is the
distance between the point charges, eps0 is the electric constant (vacuum permittivity). Molecular units are
assumed.
"""

from potentials.potential import Potential
from util import units_mu
import math

class CoulombPotential(Potential):
    """Energy and force due to interacting point charges."""

    def __init__(self, q: float):
        """Constructor
        :param q: Charge value."""
        self.q = q

    def u(self, r: float, param: dict) -> float:
        q = float(param['q'])
        return self.q * q / (4.0 * math.pi * units_mu.E0 * r)

    def du_dr(self, r: float, param: dict) -> float:
        q = float(param['q'])
        return -self.q * q / (4.0 * math.pi * units_mu.E0 * r * r)

