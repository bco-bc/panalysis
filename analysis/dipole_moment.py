"""Dipole moment M(t) fluctuations: dipole moment autocorrelation function (ACF),
dipole moment orientation probability density function (orientation relative to positive z-axis),
average dipole moment (components, strength (norm)), group dipole moment profile along a direction
(particular component of group dipole moment).
"""
import math

import util
from analysis.analyzer import Analyzer
from particle.particle_system import ParticleSystem
import numpy as np
from queue import Queue
import logging


class DipoleMoment(Analyzer):

    def __init__(self,
                 dt: float,
                 t_max: float,
                 pbc: [],
                 bin_size,
                 box: np.ndarray,
                 direction: str = 'z',
                 direction_2: str = 'x',
                 group_bin_size: float = 0.1):
        """Constructor
        :param dt Time difference between states
        :param t_max Length of time interval
        :param pbc Directions along which PBC should be applied. Subsets of {0,1,2}
        :param bin_size Bin size for probability density function of theta and cos(theta), where
        theta is the angle between the total dipole moment and the positive z-axis
        :param box Simulation box dimensions
        :param direction The direction along which particle group dipole moment profile should be computed.
        One of {x,y,z}. Default is 'z', so that the profile of the z-component will be determined
        :param direction_2 The particle group dipole moment profile is computed perpendicular to this
        direction. One of {x, y, z}. Default is 'x'.
        """
        self.dt = dt
        self.t_max = t_max
        self.pbc = pbc
        self.bin_size = bin_size
        self.box = box

        # Set up.
        self.counter = 0

        # ACF
        self.n_bins = int(self.t_max / self.dt)  # Number of bins for ACF.
        self.acf_bin_size = self.t_max / self.n_bins  # Adjusted bin size for ACF.
        self.m_0 = np.zeros(shape=3)  # Dipole moment at the beginning of the time interval.
        self.m_queue = Queue()  # Dipole moments throughout time interval.
        self.m_acf = np.zeros(shape=self.n_bins)  # Dipole moment ACF.
        self.m = np.zeros(shape=0)  # Dipole moment at time t.

        # Theta in [0, pi], cos(theta) in [-1,+1]
        self.n_bins_theta = int((math.pi - 2.0 * self.bin_size) / self.bin_size)
        self.n_bins_cos_theta = int((2.0 - 2.0 * self.bin_size) / self.bin_size)
        self.thetas = []
        self.cos_thetas = []

        # Dipole moment profile
        if direction == 'x':
            self.coordinate = 0
        elif direction == 'y':
            self.coordinate = 1
        else:
            self.coordinate = 2
        if direction_2 == 'x':
            self.coordinate_2 = 0
        elif direction_2 == 'y':
            self.coordinate_2 = 1
        else:
            self.coordinate_2 = 2
        self.n_bins_group_dipoles = int(self.box[self.coordinate_2] / group_bin_size)
        self.bin_size_group_dipole_2 = self.box[self.coordinate_2] / self.n_bins_group_dipoles
        self.group_dipoles = np.zeros(shape=self.n_bins_group_dipoles)
        self.number_of_groups = 0

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        t = self.counter * self.dt

        if self.counter == 1:
            self.number_of_groups = len(particle_system.groups)
            self.m_0 = particle_system.dipole_moment(self.pbc)
            q = particle_system.charge()
            logging.info(f'Total charge particle system: {q} ')

        # Total dipole moment.
        m_now = particle_system.dipole_moment(self.pbc)
        # print(f'counter {self.counter}: Dipole moment: {m_now} and norm {np.linalg.norm(m_now)} at time {t}.')
        self.m = np.append(self.m, m_now)
        self.m_queue.put(m_now)
        if t > self.t_max:
            self.m_0 = self.m_queue.get()

        # ACF
        for n in np.arange(0, len(self.m_queue.queue)):
            m_n = self.m_queue.get()
            inner = np.inner(self.m_0, m_n)
            self.m_acf[n] += inner
            self.m_queue.put(m_n)

        # Theta
        cos_theta = util.cvector.cos_polar_angle(m_now)
        theta = math.acos(cos_theta)
        self.cos_thetas.append(cos_theta)
        self.thetas.append(theta)

        # Group dipole moment profile.
        for g in particle_system.groups:
            dipole = g.dipole_moment(particle_system.box, self.pbc)
            r = util.box.pbc_place_inside(self.box, g.position(), pbc = [0,1,2])
            index = int(r[self.coordinate_2] / self.bin_size_group_dipole_2)
            if 0 <= index < self.n_bins_group_dipoles:
                d = dipole[self.coordinate]
                self.group_dipoles[index] += d
            else:
                logging.warning('Group position coordinate outside allowed range.')
                logging.warning(f'Particle position coordinate: {r}')
                logging.warning(f'Index value: {index}')

    def results(self):
        """Returns [t, M(t)] over the full trajectory, <M(t)> averaged over the full trajectory,
        [theta, p(theta)], where theta is the angle with the positive z-axis and p is the
        probability density function for theta, [cos(theta), p(cos(theta))],
        where cos(theta) is the cos of theta and p(cos(theta)) is its probability density function,
        [t, <M(0)M(t)>/<M(0)M(0)>] is over the time interval [0,t_max], [r, m_profile], where r is the
        location along the requested direction and m_profile is the dipole moment profile in the
        requested direction, and counter, the number of states. M(t) is the total dipole moment at time t.
        """
        # ACF
        self.m_acf /= self.counter
        m_acf_0 = self.m_acf[0]
        self.m_acf /= m_acf_0
        t_m_acf = np.zeros(shape=self.m_acf.shape)
        for k in np.arange(0, self.n_bins):
            t_m_acf[k] = k * self.acf_bin_size

        # Dipole moment
        a = self.m.ravel()
        self.m = a.reshape(self.counter, 3)
        t_m = np.zeros(shape=self.counter)
        for k in np.arange(0, self.counter):
            t_m[k] = (k + 1) * self.dt
        total_m = sum(self.m[:])
        ave_m = total_m / self.counter

        # Probability density function of theta
        p_theta, p_theta_bin_edges = np.histogram(a=self.thetas,
                                                  bins=self.n_bins_theta,
                                                  range=(self.bin_size, math.pi - self.bin_size),
                                                  density=True)
        theta = np.zeros(shape=len(p_theta))
        for k in np.arange(0, len(p_theta)):
            theta[k] = p_theta_bin_edges[k]

        # Probability density function of cos(theta).
        p_cos_theta, p_cos_theta_bin_edges = np.histogram(a=self.cos_thetas,
                                                          bins=self.n_bins_cos_theta,
                                                          range=(-1.0 + self.bin_size, 1.0 - self.bin_size),
                                                          density=True)
        cos_theta = np.zeros(shape=len(p_cos_theta))
        for k in np.arange(0, len(p_cos_theta)):
            cos_theta[k] = p_cos_theta_bin_edges[k]

        # Group dipole moment profile
        m_profile = self.group_dipoles / self.counter
        r = np.arange(0, self.box[self.coordinate_2], self.bin_size_group_dipole_2)

        return ([t_m, self.m], ave_m, self.thetas, self.cos_thetas,
                [theta, p_theta], [cos_theta, p_cos_theta], [t_m_acf, self.m_acf], [r, m_profile],self.counter)
