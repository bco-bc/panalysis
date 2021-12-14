"""Dipole moment M(t) fluctuations
"""

from analysis.analyzer import Analyzer
from particle.particle_system import ParticleSystem
import numpy as np
import logging
import util


class DipoleMoment(Analyzer):

    z_axis = np.array([0.0, 0.0, 1.0])
    """Unit vector along the positive z-axis
    """

    def __init__(self, dt: float, t_max: float, cos_a_bin_size=0.01):
        """Constructor
        :param dt Time difference between states
        :param t_max Length of time interval
        :param cos_a_bin_size Bin size for probability density function of cos(alpha), where
        alpha is the angle between the total dipole moment and the positive z-axis.
        """
        self.dt = dt
        self.t_max = t_max
        self.cos_a_bin_size = cos_a_bin_size

        self.counter = 0
        self.n_bins = int(self.t_max / self.dt)
        self.bin_size = self.t_max / self.n_bins
        self.m = np.zeros(shape=0)

        # From -1 to +1
        self.cos_a_histogram = np.zeros(shape=int(2.0/cos_a_bin_size))

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1

        if self.counter == 1:
            q = particle_system.charge()
            logging.info(f'Total charge particle system: {q} ')

        # Total dipole moment.
        origin = np.array([0.0, 0.0, 0.0])
        m_total = np.zeros(shape=(1, 3))
        for p in particle_system.all:
            delta_r = util.box.pbc_distance(particle_system.box, p.r, origin)
            m_total += p.charge() * delta_r
        self.m = np.append(self.m, m_total)

        # Probability density function of cos(alpha))
        cos_a = DipoleMoment.cos_angle__(m_total)
        index = int((cos_a + 1.0) / self.cos_a_bin_size)
        self.cos_a_histogram[index] += 1

    def results(self):
        """Returns t, M(t) over the full trajectory, <M(t)>, P(cos(alpha)), where alpha is the
        angle with the positive z-axis and P is a probability density,
        as well as <M(0)M(s)>/<M(0)M(0> over the time interval [0,t_max]. M(t) is the total
        dipole moment at time t.
        """
        a = self.m.ravel()
        b = a.reshape(self.counter, 3)
        t = np.zeros(shape=self.counter)
        for k in np.arange(0, self.counter):
            t[k] = (k + 1) * self.dt
        ave_m = sum(b[:]) / self.counter

        # Average. Compute probability density function for cos(alpha)).
        p = self.cos_a_histogram / self.counter
        total = sum(p[:])
        p_normalized = p / (total * self.cos_a_bin_size)
        return t, b, ave_m, p_normalized

    @staticmethod
    def cos_angle__(m: np.ndarray) -> float:
        """Returns cos(alpha) where alpha is the angle of m with the positive z-axis.
        """
        return  np.inner(DipoleMoment.z_axis, m) / np.linalg.norm(m)
