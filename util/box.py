"""Simulation box utilities
"""

import numpy as np
import math


def pbc_place_inside(box: np.ndarray, r_out: np.ndarray, pbc: []) -> np.ndarray:
    """Places point inside simulation box applying periodic boundary condition.
    :param box Simulation box.
    :param r_out Point possibly outside box.
    :param pbc Directions specification. Subsets of {0,1,2}.
    :return r_in Point inside simulation.
    """
    r_in = np.array([r_out[0], r_out[1], r_out[2]])
    for k in pbc:
        r_k = r_in[k]
        box_k = box[k]
        n = np.rint(math.floor(r_k/box_k))
        r_in[k] -= n * box_k
    return r_in


def pbc_distance(box: np.ndarray, r_i: np.ndarray, r_j: np.ndarray, pbc: []) -> np.ndarray:
    """Returns distance rij =|r_i - r_j| between two points applying periodic boundary condition.
    :param box Simulation box.
    :param r_i Point i.
    :param r_j Point j.
    :param pbc Directions specification. Subsets of {0,1,2}.
    :return rij = r_i - r_j
    """
    r_ij = r_i - r_j
    for k in pbc:
        box_k = box[k]
        dr = r_i[k] - r_j[k]
        ratio = dr / box_k
        n = np.rint(ratio)
        r_ij[k] = dr - n * box_k
    return r_ij
