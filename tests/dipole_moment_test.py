import numpy as np

from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory
from analysis.dipole_moment import DipoleMoment
from analysis import analyze
import matplotlib.pyplot as plt
import logging

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    fn = '/wrk3/simulation/one_water_dpd/water-specs.dat'
    catalog = read_particle_spec_catalog(fn)

    fn = '/wrk3/simploce/pt-cgmd/particles/resources/one-water.ps'
    particle_system = read_particle_sys(fn, catalog)
    fn = '/wrk3/simulation/one_water_dpd/trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)

    analyzer = DipoleMoment(dt=0.01, t_max=20, pbc=[0,1,2], bin_size=0.01)
    analyze.perform(analyzer, particle_system, trajectory, n_skip=1000)
    trajectory.close()

    t_m, ave_m, thetas, cos_thetas, theta_density, cos_theta_density, t_m_acf, n_states = analyzer.results()
    v, p_v = theta_density
    plt.plot(v, p_v, color='blue')

    plt.show()




