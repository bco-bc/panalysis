"""Interaction potential interface
"""


class Potential:
    """Defines an interface for potentials."""

    def u(self, r: float, param: dict) -> float:
        """Interaction energy U(r)
        :param r: Distance
        :param param: Parameters
        :returns Value of potential at r.
        """
        pass

    def du_dr(self, r: float, param: dict) -> float:
        """Returns dU(r)/dr. Note that force F = -dU(r)/dr.
        :param r: Distance
        :param param: Parameters
        :returns Value of derivative of potential at r.
        """
        pass
