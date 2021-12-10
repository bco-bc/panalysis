import numpy as np
from particle.particle_spec_catalog import read_particle_spec_catalog
from particle.particle_system import read_particle_sys
import simulation.trajectory
from analysis.msd import MSD
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import logging


if __name__ == '__main__':
    fn = '/localdisk/resources/particles-specs.dat'
    catalog = read_particle_spec_catalog(fn)

    fn = '/home/juffer/simulations/electrolyte/electrolyte.ps'
    particle_system = read_particle_sys(fn, catalog)
    fn = '/home/juffer/simulations/electrolyte/trajectory.dat'
    trajectory = simulation.trajectory.Trajectory(fn)

    analyzer = MSD(dt=20.0, t_max=500.0, spec=catalog.find('Na+'))
    while trajectory.next(particle_system):
        analyzer.perform(particle_system)
    trajectory.close()

    t, msd = analyzer.results()
    x = t.reshape((-1, 1))
    model = LinearRegression()
    model.fit(x, msd)
    r_sq = model.score(x, msd)
    print("Coefficients of determination:", r_sq)
    a = model.coef_[0]
    b = model.intercept_
    print("Intercept: ", b)
    print("Slope: ", a)
    prediction = np.zeros(shape=len(t))
    index = 0
    for v in t:
        prediction[index] = a * v + b
        index += 1

    plt.plot(t, msd, color='red')
    plt.plot(t, prediction, color='blue')
    plt.show()
