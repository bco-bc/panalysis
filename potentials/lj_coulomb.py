"""Several implementations for calculating the Coulomb and/or Lennard Jones
potentials between a pair of particles"""

import numpy as np
import math
from scipy import special

e0 = 0.000572766  # 'Molecular units', units of (mol e^2)/(kJ nm).
pi = np.pi


def lj_coulomb(r, param):
    """Computes the Lennard-Jones (LJ) plus the Coulomb potentials at the given distance and potentials parameters.

    Parameters
    ----------
        r : float or ndarray.
            Distance or distances (nm).
        param : tuple
            Must hold 'q1', 'q2', 'eps', 'C12', and 'C6', in that order. All float values. The latter two are the
            C12 (kJ nm^12/mol) and C6 (kJ nm^6/mol) LJ potentials parameters. The first three are the two charge
             values (e) and the relative permittivity (dielectric constant) of the Coulomb (electrostatic)
             potentials.

    Returns
    -------
        tuple
            Total potentials energy, electrostatic part and LJ part. All in kJ/mol.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    c12 = param[3]
    c6 = param[4]
    el = q1 * q2 / (4.0 * pi * e0 * eps * r)
    r2 = r * r
    r6 = r2 * r2 * r2
    r12 = r6 * r6
    lj = c12 / r12 - c6 / r6
    total = lj + el
    return total, el, lj


def coulomb_cutoff(r, param):
    """Returns the Coulomb potentials, zero beyond the cutoff distance. This implementation is referred to as Coulomb.

    Parameters
    ----------
        r : ndarray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', and 'rc', that is two charge values, the relative permittivity
            (dielectric constant), and a cutoff distance (nm).

    Returns
    -------
        ndarray of floats
            Electrostatic potentials energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    el = np.zeros(shape=0)
    for rij in r:
        if rij <= rc:
            v = q1 * q2 / rij
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def shifted_force(r, param):
    """Returns the shifted force ("SF") electrostatic energy by Levitt et al., Comput. Phys. Commun.
    1995, 91, 215−231. See also Kale, S.; Herzfeld, J. J. Chem. Theory Comput. 2011, 7, 3620−3624.

    Parameters
    ----------
        r : ndarray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', and 'rc', that is two charge values, the relative permittivity
            (dielectric constant), and a cutoff distance (nm).

    Returns
    -------
        ndarray
            Electrostatic potentials energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    el = np.zeros(shape=0)
    rc2 = rc * rc
    for rij in r:
        if rij <= rc:
            # Levitt et al, Eq (6) (top line, C = q1*q2).
            v = q1 * q2 * (1.0 / rij - 1.0 / rc + (rij - rc) / rc2)
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def damped_shifted_force(r, param):
    """Returns the damped shifted force ("DSF") electrostatic energy by Zahn et al., J. Phys. Chem. B 2002, 106, 10725.
    See also Fennell et al., J. Chem. Phys. v. 124, p. 234104, 2006.

    Parameters
    ----------
        r : ndarray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', 'rc', and 'alpha', that is two charge values,
            the relative permittivity (dielectric constant), cutoff distance (nm),
            and a damping parameter (nm^-1).

    Returns
    -------
        ndarray of floats
            Electrostatic potentials energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    alpha = param[4]
    rc2 = rc * rc
    alpha2 = alpha * alpha
    erfc_a_rc = special.erfc(alpha * rc)
    f2 = erfc_a_rc / (rc * rc)
    f3 = 2.0 * alpha * math.exp(-alpha2 * rc2) / (math.sqrt(pi) * rc)
    f = f2 + f3
    el = np.zeros(shape=0)
    for rij in r:
        if rij <= rc:
            # From Eq (13) in Zahn et al and Eq (5) in  Fennell et al.
            v = q1 * q2 * (special.erfc(alpha * rij) / rij - f * (rij - rc))
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def shifted_force_gradient(r, param):
    """Returns the shifted force gradient ("SFG") electrostatic energy by
    Kale, S.; Herzfeld, J. J. Chem. Theory Comput. 2011, 7, 3620−3624.

    Parameters
    ----------
        r : ndarray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps', and 'rc', that is two charge values, the relative permittivity
            (dielectric constant), and a cutoff distance (nm).

    Returns
    -------
        ndarray of floats
            Electrostatic potentials energies at given distances.
    """
    q1 = param[0]
    q2 = param[1]
    eps = param[2]
    rc = param[3]
    rc2 = rc * rc
    rc3 = rc2 * rc
    # rc4 = rc3 * rc
    el = np.zeros(shape=0)
    for rij in r:
        if rij <= rc:
            rij_rc = rij - rc
            rij_rc_2 = rij_rc * rij_rc
            # rij_rc_3 = rij_rc_2 * rij_rc
            # Eq (5) in Kale et al.
            v = q1 * q2 * (1.0 / rij - 1.0 / rc + rij_rc / rc2 - rij_rc_2 / rc3)
            # v = q1 * q2 * (1.0 / rij - 1.0 / rc + rij_rc / rc2 - rij_rc_2 / rc3 + rij_rc_3 / rc4)
            v /= (4.0 * pi * e0 * eps)
            el = np.append(el, v)
        else:
            el = np.append(el, 0.0)
    return el


def reaction_field(r, param) -> (np.ndarray, np.ndarray):
    """Calculates the electrostatic potentials based on the reaction field (RF) approach as listed in
    Riniker and van Gunsteren, J. Chem. Phys. 134, 084110, 2011.

    Parameters
    ----------
        r : ndarray of floats
            Distances (nm)
        param : tuple
            Must hold 'q1', 'q2', 'eps_cs', 'eps_rf, and 'rc', that is two charge values, the relative permittivity
            within the cutoff distance and the relative permittivity outside the cutoff distance, and a cutoff distance
            (nm).

    Returns
    -------
        two ndarray of floats
            Electrostatic potentials energies and associated force at given distances.
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
    c_rf = ((2.0 * eps_cs - 2.0 * eps_rf) * f1 - f2) / ((eps_cs + 2.0 * eps_rf) * f1 + f2)
    print('Reaction field: C_rf = ', c_rf)
    print('Reaction field: eps_cs = ', eps_cs)
    print('Reaction field: eps_rf = ', eps_rf)
    el = np.zeros(0)
    forces = np.zeros(0)
    factor = 4.0 * pi * e0 * eps_cs
    c1 = q1 * q2 / factor
    for rij in r:
        if rij <= rc:
            rij2 = rij * rij
            c = c1 / rij
            dc_dr = -c1 / rij2
            rf = -c1 * (0.5 * c_rf * rij2 / rc3 + (1.0 - 0.5 * c_rf) / rc)
            drf_dr = -c1 * c_rf * rij / rc3
            v = c + rf
            force = -(dc_dr + drf_dr)
            el = np.append(el, v)
            forces = np.append(forces, force)
        else:
            el = np.append(el, 0.0)
            forces = np.append(forces, 0.0)
    return el, forces


def gaussian_density(r, params) -> (np.ndarray, np.ndarray):
    """Electric potential due to a Gaussian charge density of the form
    rho(r)=(Q / pi^3/2 * sigma) * exp(-r^2/sigma^2), where sigma (width of the density) and Q (total charge)
    are given.
    param r: Distances relative to charge density for which potential must be computed
    :param params A tuple with parameters sigma and Q
    """
    sigma = params[0]
    q = params[1]
    gauss = np.zeros(0)
    coulomb = np.zeros(0)
    sqrt_pi = math.sqrt(math.pi)
    f = q / (4.0 * math.pi * e0)
    for rij in r:
        if rij < 1.0e-03:
            c = 100.0 * f
            g = 2.0 / (sqrt_pi * sigma) * f
        else:
            c = f / rij
            g = math.erf(rij / sigma) * c
        gauss = np.append(gauss, g)
        coulomb = np.append(coulomb, c)
    return coulomb, gauss


def solid_sphere_density(r, params) -> (np.ndarray, np.ndarray):
    """Electric potential due to a uniformly charged solid sphere, that is rho(r)=rho for r<=R and
    rho(r)=0 for r>R, where R is the radius of the sphere.
    :param r: Distances relative to charge density for which potential must be computed
    :param params: A tuple with parameters R and Q, where Q is the total charge of the sphere
    :return Coulomb and solid sphere potential
    """
    radius = params[0]
    q = params[1]
    solid = np.zeros(0)
    coulomb = np.zeros(0)
    f = q / (4.0 * math.pi * e0)
    for rij in r:
        if rij < 1.0e-03:
            c = 100.0 * f
        else:
            c = f / rij
        if rij <= radius:
            s1 = 3.0 - (rij * rij) / (radius * radius)
            s2 = f / (2.0 * radius)
            s = s2 * s1
        else:
            s = f / rij
        solid = np.append(solid, s)
        coulomb = np.append(coulomb, c)
    return coulomb, solid
