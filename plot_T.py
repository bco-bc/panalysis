"""Extracts the temperature from a simulation data output file, and generates a plot versus time, along with
requested reference temperature and the time average.
"""

import util
import sys
import numpy as np
import matplotlib.pyplot as plt
import logging


def usage():
    print()
    print('Usage: python plot_T.py [OPTION [value]]...')
    print('Description: Plots temperature from simulation data versus a reference temperature.')
    print()
    print('Optional arguments:')
    print('-dt VALUE: VALUE is the time difference (in ps) between states in trajectory. Default is 0.002.')
    print('-fn-data FN: FN is the simulation data file. Default is \'simulation.dat\'')
    print('-T VALUE: Value is the reference temperature. Default is 298.15.')
    print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    conf = util.parse_argv.parse(sys.argv)
    print(f'Received arguments: {conf}')
    fnData = 'simulation.dat'
    try:
        fnData = conf['fn-data']
        temperature = float(conf['temperature'])
        dt = conf['dt']
    except KeyError:
        usage()
        sys.exit("Missing input arguments")

    print('Reference temperature: ', temperature)

    # Read data.
    data = np.loadtxt(fnData, usecols=(1, 8))
    t = data[:, 0]
    T = data[:, 1]
    T_ref = np.full(shape=t.shape, fill_value=temperature)
    ave = np.arange(0)
    counter = 0
    total = 0.0
    for temp in T:
        counter += 1
        total += temp
        v = total / counter
        ave = np.append(ave, v)

    # Plot
    a = ave[ave.size - 1]
    plt.plot(t, T, label=r'$T(t)$')
    label = r'$<T>$ = ' + str(a)
    plt.plot(t, ave, label=label)
    label = r'$T_{ref}$ = ' + str(temperature)
    plt.plot(t, T_ref, label=label)
    plt.xlabel('t')
    plt.ylabel('T')
    plt.legend()
    plt.show()
