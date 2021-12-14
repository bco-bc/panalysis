"""Some utilities for cartesian vectors.
"""

import numpy as np
import math


def inner(x_1: np.ndarray, x_2: np.ndarray) -> float:
    """Returns inner product between two cartesian vectors
    :param x_1 Cartesion vector with shape=(3,1) or (1,3)
    :param x_2 Cartesion vector with shape=(3,1) or (1,3)
    """
    return np.inner(x_1, x_2)


def norm(x: np.ndarray) -> float:
    """Returns norm or length of a cartesion vector
    :param x Vector with shape=(3,1) or (1,3)
    """
    v = inner(x, x)
    return math.sqrt(v)
