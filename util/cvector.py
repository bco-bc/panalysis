"""Some utilities for cartesian vectors.
"""
import math

import numpy as np


def cos_polar_angle(v: np.ndarray) -> float:
    """Returns cos(theta) where theta is the angle of the Cartesian vector v with the
    positive z-axis.
    :param v Cartesian vector with nonzero norm
    :returns cos(theta) between -1 and +1.
    """
    return v[2] / np.linalg.norm(v)


def polar_angle(v: np.ndarray) -> float:
    """Returns theta where theta is the angle of the Cartesian vector v with the
    positive z-axis.
    :param v Cartesian vector.
    :returns Theta between 0 and pi.
    """
    cos_a = cos_polar_angle(v)
    a = math.acos(cos_a)
    return a


def azimuthal_angle(v: np.ndarray) -> float:
    """Returns azimuthal angle of given Cartesian vector
    :param v Cartesian vector
    :returns Azimuthal angle between 0 and 2*pi
    """
    x = v[0]
    y = v[1]
    if x == 0.0 and y == 0.0:
        raise Exception(f'{v}: azimuthal angle not defined for this Cartesian vector')

    if x == 0.0 and y > 0.0:
        a = 0.5 * math.pi
    elif x == 0.0 and y < 0.0:
        a = -0.5 * math.pi
    else:
        a = math.atan2(y, x)
    a += math.pi
    return a
