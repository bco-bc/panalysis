"""Displacement
"""

from analysis.analyzer import Analyzer
from particle.particle_spec import ParticleSpec
from particle.particle_system import ParticleSystem
import numpy as np
from queue import Queue


class MSD(Analyzer):
    """Calculates the average mean square displacement (MSD) of particles over a given
    time interval.
    """

    def __init__(self, dt: float, t_max: float, spec: ParticleSpec):
        """Constructor
        :param dt Length of time interval between successive states
        :param t_max Length of time interval
        :param spec Particle specification identifying particle for which the displacement is computed
        """
        self.dt = dt
        self.t_max = t_max
        self.spec = spec

        self.counter = 0
        self.n_bins = int(self.t_max / self.dt)
        self.bin_size = self.t_max / self.n_bins
        self.n_specs = 0
        self.r_i = np.zeros(shape=0)                 # Positions at t in [t, t+t_max]
        self.r_queue = Queue()                       # Positions in [t, t+t_max]
        self.msd = np.zeros(shape=self.n_bins)       # MSD

    def perform(self, particle_system: ParticleSystem):
        self.counter += 1
        t = self.counter * self.dt

        if self.counter == 1:
            self.n_specs = Analyzer.number_of_specs(particle_system, self.spec)
            if self.n_specs == 0:
                raise ValueError(f'{self.spec.name}: No such particle with this specification.')
            self.r_i = Analyzer.positions(particle_system)

        r_now = Analyzer.positions(particle_system)
        self.r_queue.put(r_now)
        if t > self.t_max:
            self.r_i = self.r_queue.get()

        for n in np.arange(0, len(self.r_queue.queue)):
            r_n = self.r_queue.get()
            for k in np.arange(0, len(r_n)):
                particle = particle_system.all[k]
                if particle.spec.name == self.spec.name:
                    delta_r = r_n[k] - self.r_i[k]
                    self.msd[n] += np.inner(delta_r, delta_r)
            self.r_queue.put(r_n)

    def results(self) -> (np.ndarray, np.ndarray):
        """Returns t and msd(t)
        """
        self.msd /= (self.n_specs * self.counter)
        t = np.zeros(shape=self.msd.size)
        for k in np.arange(0, self.msd.size):
            t[k] = k * self.bin_size
        return t, self.msd
