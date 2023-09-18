"""Simple utility to parse program arguments
"""

import sys
import logging


def parse(argv: sys.argv) -> dict:
    """Parses program options. Read options and associated values. E.g.,
    '-s particle-specs.dat', where the option is '-s' and the associated value is
    'particle-specs.dat'
    """

    # Defaults
    configuration = {'fn-trajectory': 'trajectory.dat',
                     'fn-data': 'simulation.dat',
                     'fn-particle-specs': 'particle-specs.dat',
                     'temperature': 298.15,
                     'distance' : 0.0,
                     'E0': 0.0,
                     'direction': 'z',
                     'second-direction': 'x',
                     'pbc': [0, 1, 2],
                     'bin-size': 0.1,
                     'exclude': None,
                     'dt': 0.002,
                     'n-skip' : 0,
                     'id-1': '',
                     'id-2': '',
                     'spec-1': '',
                     'spec-2': '',
                     't-max': 0.0}

    # Default settings.
    counter = 1
    while counter < len(argv):
        item = argv[counter]
        if item == '-s':
            # File name particle specifications.
            counter += 1
            configuration['fn-particle-specs'] = argv[counter]
            entry = configuration['fn-particle-specs']
            logging.info(f'Particle specifications file name: {entry}')
        if item == '-tr':
            # File name trajectory.
            counter += 1
            configuration['fn-trajectory'] = argv[counter]
            entry = configuration['fn-trajectory']
            logging.info(f'Trajectory file name: {entry}')
        if item == '-fn-data':
            # File name simulation data
            counter += 1
            configuration['fn-data'] = argv[counter]
            entry = configuration['fn-data']
            logging.info(f'Simulation data file name: {entry}')
        if item == '-ps':
            # File name particle system.
            counter += 1
            configuration['fn-particle-system'] = argv[counter]
            entry = configuration['fn-particle-system']
            logging.info(f'Particle system file name: {entry}')
        if item == '-dt':
            # Time difference between entries in trajectory.
            counter += 1
            configuration['dt'] = float(argv[counter])
            entry = configuration['dt']
            logging.info(f'Time difference between trajectory entries: {entry}')
        if item == '-t-max':
            # Length of time interval.
            counter += 1
            configuration['t-max'] = float(argv[counter])
            entry = configuration['t-max']
            logging.info(f'Length time interval: {entry}')
        if item == '-r-max':
            # Length of distance interval.
            counter += 1
            configuration['r-max'] = float(argv[counter])
            entry = configuration['r-max']
            logging.info(f'Length distance interval: {entry}')
        if item == '-E0':
            # External electric field strength
            counter += 1
            configuration['E0'] = float(argv[counter])
            entry = configuration['E0']
            logging.info(f"Strength external electric field: {entry}")
        if item == '-d':
            # Direction of something.
            counter += 1
            configuration['direction'] = argv[counter]
            entry = configuration['direction']
            logging.info(f"Direction: {entry}")
        if item == '-d2':
            # Direction of something.
            counter += 1
            configuration['other-direction'] = argv[counter]
            entry = configuration['other-direction']
            logging.info(f"Second direction: {entry}")
        if item == '-T':
            # Temperature
            counter += 1
            configuration['temperature'] = float(argv[counter])
            entry = configuration['temperature']
            logging.info(f'Temperature: {entry}.')
        if item == '-R':
            # Distance
            counter += 1
            configuration['distance'] = float(argv[counter])
            entry = configuration['distance']
            logging.info(f'Distance: {entry}.')
        if item == '-spec':
            # Particle specification name.
            counter += 1
            configuration['spec'] = argv[counter]
            entry = configuration['spec']
            logging.info(f'Particle specification name: {entry}.')
        if item == '-spec-1':
            # Particle specification name.
            counter += 1
            configuration['spec-1'] = argv[counter]
            entry = configuration['spec-1']
            logging.info(f'Particle specification name #1: {entry}.')
        if item == '-spec-2':
            # Particle specification name.
            counter += 1
            configuration['spec-2'] = argv[counter]
            entry = configuration['spec-2']
            logging.info(f'Particle specification name #1: {entry}.')
        if item == '-id-1':
            # Particle identifier
            counter += 1
            configuration['id-1'] = argv[counter]
            entry = configuration['id-1']
            logging.info(f'Particle identifier #1: {entry}.')
        if item == '-id-2':
            # Particle identifier
            counter += 1
            configuration['id-2'] = argv[counter]
            entry = configuration['id-2']
            logging.info(f'Particle identifier #2: {entry}.')
        if item == '-pbc-1':
            # PBD in 1D only.
            counter += 1
            direction = argv[counter]
            if direction == 'x':
                configuration['pbc'] = [0]
            elif direction == 'y':
                configuration['pbc'] = [1]
            else:
                configuration['pbc'] = [2]
        if item == '-bin-size':
            # Bin size for some histogram.
            counter += 1
            configuration['bin-size'] = argv[counter]
            entry = configuration['bin-size']
            logging.info(f'Bin size histogram: {entry}')
        if item == '--dielectric-constant':
            counter += 1
            configuration['epsilon'] = argv[counter]
            entry = configuration['epsilon']
            logging.info(f'Dielectric constant: {entry}')
        if item == '--radius':
            counter += 1
            configuration['radius'] = argv[counter]
            entry = configuration['radius']
            logging.info(f'Radius: {entry}')
        if item == '-exclude':
            counter += 1
            configuration['exclude'] = argv[counter]
            entry = configuration['exclude']
            logging.info(f'Excluding {entry}')
        if item == '-n-skip':
            counter += 1
            configuration['n-skip'] = argv[counter]
            entry = configuration['n-skip']
            logging.info(f'Skipping the first {entry} state(s) from trajectory')
        counter += 1

    return configuration
