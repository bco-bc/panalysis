import numpy as np

e0 = 0.000572766  # 'Molecular units', units of (mol e^2)/(kJ nm).
pi = np.pi


def potential(r, param):
    """Computes the Lennard-Jones (LJ) plus the Coulomb interaction at the given distance and interaction parameters.

    Parameters
    ----------
        r : float
            Distance (nm).
        param : tuple
            Must hold 'q1', 'q2', 'eps', 'C12', and 'C6', in that order. All float values. The latter two are the
            C12 kJ nm^12/mol) and C6 (kJ nm^6/mol) LJ interaction parameters. The first three are the charges (e)
            and the relative permittivity (dielectric constant) of the Coulomb (electrostatic) interaction.

    Returns
    -------
        potential(r, param) : tuple
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
