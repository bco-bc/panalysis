"""Gaussian interaction potential U(r) between two overlapping charge densities.
U(r) = erf((1.0/(s1^2 * s2^2))^1/2 * r) * Q1 * Q2 / (4.0 * pi * eps0 * r) where s1 is sigma of this density and s2 is
sigma of the other density, Q1 is the total charge of this density and Q2 is the total charge of the other density,
r is the distance between the densities and eps0 is the electric constant (vacuum permittivity). Erf(x) is the error
function. Molecular units are assumed.
"""

from potentials.potential import Potential
from util import units_mu
import math


class GaussianPotential(Potential):
    """Energy and force due to overlapping Gaussian charge densities.
    """

    def __init__(self, sigma: float, q: float):
        """
        Constructor.
        :param sigma: Width of the density.
        :param q: Total charge of the density
        """
        self.sigma = sigma
        self.q = q

    def u(self, r: float, param: dict) -> float:
        """Interaction energy U(r) between two Gaussian densities.
        :param r: Distance
        :param param: Parameters. Must hold total charge ("q") and the width ("sigma") of the other Gaussian density
        :returns U(r)"""
        sigma = float(param['sigma'])
        q = float(param['q'])
        s = self.s_(sigma)
        x = s * r
        erf_x = math.erf(x)
        return erf_x * self.q * q / (4.0 * math.pi * units_mu.E0 * r)

    def du_dr(self, r: float, param: dict) -> float:
        """Derivative of interaction energy U(r) for two overlapping Gaussian densities.
        :param r: Distance
        :param param: Parameters. Must hold total charge ("q") and the width ("sigma") of the other Gaussian density
        :returns dU(r)/dr"""
        sigma = float(param['sigma'])
        q = float(param['q'])
        s = self.s_(sigma)
        x = s * r
        e = math.exp(-x * x)
        x1 = 2.0 / math.sqrt(math.pi) * e * x
        x2 = math.erf(x)
        return (x1 - x2) * self.q * q / (4.0 * math.pi * units_mu.E0 * r * r)

    def s_(self, sigma: float) -> float:
        return math.sqrt(1.0 / (self.sigma * self.sigma + sigma * sigma))
