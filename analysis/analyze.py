"""Reads states from the trajectory and passes them onto the particle system to be analyzed by
    an analyzer
"""

from simulation.trajectory import Trajectory
from analysis.analyzer import Analyzer
from particle.particle_system import ParticleSystem

def perform(analyzer: Analyzer, particle_system: ParticleSystem, trajectory: Trajectory, n_skip: int = 0):
    """Reads states from the trajectory and passes them onto the particle system to be analyzed by the
    analyzer
    :param analyzer Analyzes states assigned to the particle system
    :param particle_system States are assigned to this particle system
    :param trajectory States
    :param n_skip Number of states to skip from the trajectory before passing them on
    """
    counter = 0
    while trajectory.next(particle_system):
        counter += 1
        if counter > n_skip:
            analyzer.perform(particle_system)
