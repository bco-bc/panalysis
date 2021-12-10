"""Velocity autocorrelation function.
"""

from analysis.analyzer import Analyzer
from particle.particle_spec import ParticleSpec
from particle.particle_system import ParticleSystem
import numpy as np
from queue import Queue


class VAF(Analyzer):
    """Computes the velocity autocorrelation function (VAF) for a specific particle of a given
    time interval.
    """

    def __init__(self, dt: float, t_max: float, spec: ParticleSpec):
        """Constructor
        :param dt Time difference between states.
        :param t_max Length of time interval.
        :param spec Particle specification identifying particle for which the VAF is computed.
        """
        self.dt = dt
        self.t_max = t_max
        self.spec = spec

        self.counter = 0
        self.n_bins = int(self.t_max / self.dt)
        self.bin_size = self.t_max / self.n_bins
        self.n_specs = 0
        self.v_i = np.zeros(shape=0)               # Velocities at t in [t, t+t_max]
        self.v_queue = Queue()                     # Velocities in [t, t+t_max]
        self.vaf = np.zeros(shape=self.n_bins)     # VAF

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        t = self.counter * self.dt

        if self.counter == 1:
            self.n_specs = Analyzer.number_of_specs(particle_system, self.spec)
            self.v_i = Analyzer.velocities(particle_system)

        v_now = Analyzer.velocities(particle_system)
        self.v_queue.put(v_now)
        if t > self.t_max:
            self.v_i = self.v_queue.get()

        for n in np.arange(0, len(self.v_queue.queue)):
            v_n = self.v_queue.get()
            for k in np.arange(0, len(v_n)):
                particle = particle_system.all[k]
                if particle.spec.name == self.spec.name:
                    self.vaf[n] += np.inner(v_n[k], self.v_i[k])
            self.v_queue.put(v_n)

    def results(self) -> (np.ndarray, np.ndarray):
        """Returns t and vaf(t)/vaf(0)
        """
        self.vaf /= (self.counter * self.n_specs)
        vaf_0 = self.vaf[0]
        self.vaf /= vaf_0
        t = np.zeros(shape=self.n_bins)
        for k in np.arange(0, self.n_bins):
            t[k] = k * self.bin_size
        return t, self.vaf
