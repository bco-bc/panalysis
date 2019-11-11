"""Several implementations for calculating the Coulomb and/or Lennard Jones
interaction between a pair of particles"""

import numpy as np
import math
from scipy import special

e0 = 0.000572766  # 'Molecular units', units of (mol e^2)/(kJ nm).
pi = np.pi


def potential(r, param):
    """Computes the Lennard-Jones (LJ) plus the Coulomb interaction at the given distance and interaction parameters.

    Parameters
    ----------
        r : float or nparray.
            Distance or distances (nm).
        param : tuple
            Must hold 'q1', 'q2', 'eps', 'C12', and 'C6', in that order. All float values. The latter two are the
            C12 kJ nm^12/mol) and C6 (kJ nm^6/mol) LJ interaction parameters. The first three are the two charge
             values (e) and the relative permittivity (dielectric constant) of the Coulomb (electrostatic)
             interaction.

    Returns
    -------
        tuple
            Total interaction energy, electrostatic part and LJ part. All in kJ/mol.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    C12 = param[3]
    C6 = param[4]
    el = q1 * q2 / (4.0 * pi * e0 * eps * r)
    r2 = r * r
    r6 = r2 * r2 * r2
    r12 = r6 * r6
    lj = C12 / r12 - C6 / r6
    total = lj + el
    return total, el, lj


def coulomb_cutoff(r, param):
    """Returns the Coulomb interaction, zero beyond the cutoff distance. This implementation is referred to as RC.

    Parameters
    ----------
        r : nparray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', and 'rc', that is two charge values, the relative permittivity
            (dielectric constant), and a cutoff distance (nm).

    Returns
    -------
        nparray of floats
            Electrostatic interaction energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    el = np.arange(0)
    for rij in r:
        if rij <= rc:
            v = q1 * q2 / rij
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def shifted_force(r, param):
    """Returns the shifted force (SF) electrostatic energy by Levitt et al., Comput. Phys. Commun. 1995, 91, 215−231.

    Parameters
    ----------
        r : nparray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', and 'rc', that is two charge values, the relative permittivity
            (dielectric constant), and a cutoff distance (nm).

    Returns
    -------
        nparray
            Electrostatic interaction energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    el = np.arange(0)
    rc2 = rc * rc
    for rij in r:
        if rij <= rc:
            v = q1 * q2 * (1.0 / rij - 1.0 / rc + (rij - rc) / rc2)
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def damped_shifted_force(r, param):
    """Returns the damped shifted force (DSF) electrostatic energy by Zahn et al, J. Phys. Chem. B 2002, 106, 10725.

    Parameters
    ----------
        r : nparray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', 'rc', and 'alpha', that is two charge values, the relative permittivity
            (dielectric constant), cutoff distance (nm), and a damping parameter (nm^-1).

    Returns
    -------
        nparray of floats
            Electrostatic interaction energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    alpha = param[4]
    rc2 = rc * rc
    alpha2 = alpha * alpha
    erfc_a_rc = special.erfc(alpha * rc)
    f1 = erfc_a_rc / rc
    f2 = erfc_a_rc / (rc * rc)
    f3 = 2.0 * alpha * math.exp(-alpha2 * rc2) / (math.sqrt(pi) * rc)
    f = f2 + f3
    el = np.arange(0)
    for rij in r:
        if rij <= rc:
            v = q1 * q2 * (special.erfc(alpha * rij) / rij - f1 + f * (rij - rc))
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def shifted_force_3nd_derivative(r, param):
    """Returns the shifted force electrostatic energy by Kale, S.; Herzfeld, J. J. Chem. Theory Comput. 2011, 7, 3620−
3624.

    Parameters
    ----------
        r : nparray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', and 'rc', that is two charge values, the relative permittivity
            (dielectric constant), and a cutoff distance (nm).

    Returns
    -------
        nparray of floats
            Electrostatic interaction energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    rc2 = rc * rc
    rc3 = rc2 * rc
    rc4 = rc3 * rc
    el = np.arange(0)
    for rij in r:
        if rij <= rc:
            rij_rc = rij - rc
            rij_rc_2 = rij_rc * rij_rc
            rij_rc_3 = rij_rc_2 * rij_rc
            v = q1 * q2 * (1.0 / rij - 1.0 / rc + rij_rc / rc2 - rij_rc_2 / rc3 + rij_rc_3 / rc4)
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def reaction_field(r, param):
    """Calculates the electrostatic interaction based on the reactian field (RF) approach as listed in
    Riniker and van Gunsteren, J. Chem. Phys. 134, 084110, 2011.

    Parameters
    ----------
        r : nparray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps_cs', 'eps_rf, and 'rc', that is two charge values, the relative permittivity
            within the cutoff distance and the relative permittivit outside the cutoff distance, and a cutoff distance
            (nm).

    Returns
    -------
        nparray of floats
            Electrostatic interaction energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps_cs = param[2]
    eps_rf = param[3]
    rc = param[4]
    rc2 = rc * rc
    rc3 = rc2 * rc
    kappa = param[5]
    kappa2 = kappa * kappa
    f1 = 1.0 + kappa * rc
    f2 = eps_rf * kappa2 * rc2
    C_rf = ((2.0 * eps_cs - 2.0 * eps_rf) * f1 - f2) / ((eps_cs + 2.0 * eps_rf) * f1 + f2)
    print('Reaction field: C_rf = ', C_rf)
    el = np.arange(0)
    for rij in r:
        if rij <= rc:
            rij2 = rij * rij
            v_c = q1 * q2 / rij
            v_rf = - q1 * q2 * (0.5 * C_rf * rij2 / rc3 + (1.0 - 0.5 * C_rf) / rc)
            v = v_c + v_rf
            v /= (4.0 * pi * e0 * eps_cs)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el
