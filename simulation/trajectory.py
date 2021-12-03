"""Trajectory
"""

import numpy as np
from particle.particle_system import ParticleSystem


class Trajectory():
    """Represents a trajectory produced by a simulation."""

    def __init__(self, fn: str = 'trajectory.dat'):
        """Constructor.
        fn : File name of trajectory file.
        """
        self.file = open(fn, 'r')

    def next(self, particle_system: ParticleSystem) -> bool:
        """Reads next state from simulation and assigns it to the given particle system.
        :return false if all states were read.
        """
        # Read next line from trajectory
        line = self.file.readline()
        items = line.split()

        # If empty, return.
        if len(items) == 0:
            return False

        counter = 0
        for p in particle_system.all:
            p.r = np.array([float(items[counter]), float(items[counter+1]), float(items[counter+2])])
            p.v = np.array([float(items[counter+3]), float(items[counter+4]), float(items[counter+5])])
            if p.spec.protonatable:
                pass
            counter += 6

        # Not empty.
        return True

    def close(self):
        """Close the trajectory file."""
        self.file.close()
