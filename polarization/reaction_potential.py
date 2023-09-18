"""Reaction potential inside a sphere of given radius in a dielectric medium with given dielectric constant.
"""

import math

import numpy as np

import util.units_mu as mu


def image_approximation_q(Q: float, conf: dict) -> tuple:
    """Calculates the reaction potential due to a single charge at the charge's location on the z-axis
    from the center to close to the surface inside a spherical cavity in a dielectric medium according to the
    image approximation. The reaction potential units are V.
    by Friedman, Molec. Phys., v. 29, 1533-1543, 1975.
    :param Q Charge value
    :param conf Parameters, the radius of the sphere ('radius') and the dielectric constant of medium
    outside the sphere ('eps'). The dielectric constant should be large, e.g., 80
    return z, rp (position in sphere and reaction potential at z)
    """
    a = float(conf['radius'])
    eps = float(conf['epsilon'])
    bin_size = float(conf['bin-size'])
    a2 = a * a
    factor = (1.0 - eps) * a / (eps + 1.0)
    # In
    four_pi_e0 = 4.0 * math.pi * mu.E0
    n = int(a / bin_size)
    z = np.zeros(shape=n)
    r_p = np.zeros(shape=n)
    bin_size = a / n
    for i in np.arange(0, n):
        # r is also the charge position.
        r = 0.0001 + i * bin_size
        r_image = a2 / r               # Eq 4b in Friedman, 1975.
        q_image = factor * Q / r       # Eq 4a in Friedman, 1975.
        z[i] = r
        # Potential is kJ/(mol e)
        p = q_image / (four_pi_e0 * (r_image - r))
        r_p[i] = p / mu.V_to_kJ_mol_e  # In V.
    # Return position in sphere and reaction potentials at charge location.
    return z, r_p


def image_approximation(r: np.ndarray, r_Q: np.ndarray, Q: np.ndarray, conf: dict) -> np.ndarray:
    """Calculates the reaction potential due to a set of charges at specified points inside a spherical cavity in a
    dielectric medium according to the image approximation by Friedman, Molec. Phys., v. 29, 1533-1543, 1975.
    The reaction potential units are V.
    :param r Points inside the sphere at which reaction potential is to be calculated expressed relative
    to the center of the sphere
    :param r_Q Location charges inside sphere expressed relative to the center of the sphere
    :param Q Charge values
    :param conf Parameters, the radius of the sphere ('radius') and the dielectric constant of medium
    outside the sphere ('eps'). The dielectric constant should be large, e.g., 80.
    :returns Reaction potentials at r
    """
    a = float(conf['radius'])
    eps = float(conf['epsilon'])
    a2 = a * a
    factor = (1 - eps) * a / (eps + 1)
    four_pi_e0 = 4.0 * math.pi * mu.E0
    rp = np.zeros(shape=len(r))
    for i in np.arange(0, len(r)):
        r_i = r[i]
        for j in np.arange(0, len(Q)):
            r_q = r_Q[j]
            q = Q[j]
            norm_r_q = np.linalg.norm(r_q)
            r_image = a2 * r_q / norm_r_q
            q_image = factor * q / norm_r_q
            distance = np.linalg.norm(r_i - r_image)
            # Potential is kJ/(mol e)
            p = q_image / (four_pi_e0 * distance)
            # Convert to V
            rp[i] += p / mu.V_to_kJ_mol_e
    # Return reaction potentials at r.
    return rp
