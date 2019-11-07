"""Polarizable CG water model. See Riniker and van Gunsteren, J. Chem. Phys. 134, 084110, 2011.
"""

import numpy as np
import math
from scipy import special

K_b = 2.0e+6      # Force constant of bond between CG and DP.
r_CW_DP = 0.2     # nm.
C12 = 1.298e-3    # CW - CW interaction.
C6 = 0.088        # CW - CW interaction.


def intra_dw_dp(r):
    """Calculates the bonded interaction between a DW and DP particle.

    Parameters
    ----------
        r : nparray of floats
            Holds distances between CW and CW.

    Returns
    -------
        intra_dw_dp(r) : nparray of floats.
            PInteraction nergy at the given distances.
    """
    epot = np.arange(0)
    for rij in r:
        d = r - r_CW_DP
        if d > 0.0:
            d4 = d * d * d * d;
            v = 0.5 * K_b * d4;
            epot = np.append(epot, v)
        else:
            epot = np.append(epot, 0.0)
    return epot
