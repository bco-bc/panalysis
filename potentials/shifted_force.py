"""Calculates the shifted force potentials for a given potentials.
"""
from potentials.potential import Potential
import numpy as np


def shifted_force(p: Potential, r: np.ndarray, r_C: float, param: dict):
    """Calculates SF for given potential
    :param p Original potential
    :param r Distances
    :param r_C Cutoff distance
    :param param: Parameters
    """
    u_r = np.zeros(shape=len(r))
    du_dr_r = np.zeros(shape=len(r))
    u_r_c = p.u(r_C, param)
    du_dr_r_cutoff = p.du_dr(r_C, param)
    for k in np.arange(0, len(r)):
        r_k = r[k]
        # Shift for r <= r_C
        if r_k <= r_C:
            # Original potential values
            u_r_k = p.u(r_k, param)
            du_dr_r_k = p.du_dr(r_k, param)
            # Shift
            u_r[k] = u_r_k - u_r_c - du_dr_r_cutoff * (r_k - r_C)
            du_dr_r[k] = du_dr_r_k - du_dr_r_cutoff
        else:
            u_r[k] = 0.0
            du_dr_r[k] = 0.0
    return r, u_r, du_dr_r
