"""Numerical values of physical constants in 'molecular units' (MU units). These are derived from values in SI units.
Units for MU are:
Time: ps,
Distance, position: nm (= 10^-9 m),
Velocity: nm/ps (1 ps = 10^-12 s),
Mass: u (unified atomic mass unit, 1 u = 1.66054e-27 kg)
Momentum: (u nm)/ps
Energy: kJ/mol = (u nm^2)/(ps^2)
Force: kJ/(mol nm) = (u nm)/(ps^2)
Charge: e (= 1.6021766208e-19 C).
Electric potential: kJ/(mol e)
Electric field: kJ/(mol e nm)
See also: Berendsen, H. J. C, Simulating the physical world. Cambridge University Press, 2007 (p. xv - xxvii).
"""

import scipy.constants as constants

eV = 1.602176634e-19
"""eV expressed in J."""

E0 = constants.epsilon_0 / (constants.elementary_charge * 1.0e+09) * 1.0e+03 / (eV * constants.N_A)
"""Vacuum permittivity in (mol e^2)/(kJ nm)
"""

V_to_kJ_mol_e = constants.N_A * constants.elementary_charge / 1.0e+03
"""Conversion factor for V to kJ/(mol e)
"""