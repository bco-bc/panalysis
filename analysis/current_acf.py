"""Computes the autocorrelation of the total electric current. The current is defined as
J = sum_i(j_i), where j_i = q_i * v_i, and i refers to particle i and v_i is velocity.
Consequently, the units of J and j_i are e*nm/ps or C*m/s instead of C/s=A.
"""

from analysis.analyzer import Analyzer
from particle.particle_system import ParticleSystem
import math
import numpy as np
from queue import Queue


class CurrentACF(Analyzer):
    """Calculates the total electric current autocorrelation function"""

    def __init__(self, dt: float, t_max: float):
        """Constructor
        :param dt Time difference between states.
        :param t_max Length of time interval.
        """
        self.dt = dt
        self.t_max = t_max

        self.counter = 0
        self.n_bins = int(self.t_max / self.dt)
        self.bin_size = self.t_max / self.n_bins
        self.current_0 = np.zeros(shape=3)                # Current at the beginning of the time interval.
        self.current_queue = Queue()                      # Currents throughout the time interval.
        self.current_acf = np.zeros(shape=self.n_bins)    # Current ACF.

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        t = self.counter * self.dt

        if self.counter == 1:
            self.current_0 = CurrentACF.current__(particle_system)

        current_now = CurrentACF.current__(particle_system)
        self.current_queue.put(current_now)
        if t > self.t_max:
            self.current_0 = self.current_queue.get()

        for n in np.arange(0, len(self.current_queue.queue)):
            current_n = self.current_queue.get()
            inner = np.inner(self.current_0, current_n)
            self.current_acf[n] += inner
            self.current_queue.put(current_n)

    def results(self) -> (np.ndarray, np.ndarray):
        """Returns t and J_acf(t)
        """
        self.current_acf /= self.counter
        current_acf_o = self.current_acf[0]
        self.current_acf /= current_acf_o
        t = np.zeros(shape=self.current_acf.size)
        for k in np.arange(0, self.current_acf.size):
            t[k] = k * self.bin_size
        return t, self.current_acf

    @staticmethod
    def current__(particle_system: ParticleSystem):
        total_j = np.zeros(shape=3)
        for p in particle_system.all:
            j_p = p.spec.charge * p.v
            total_j += j_p
        return total_j
