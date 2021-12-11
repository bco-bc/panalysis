"""Simple utility to parse program arguments
"""

import sys
import logging


def parse(argv: sys.argv) -> dict:
    """Parses program options. Read options and associated values, if required. E.g.,
    '-s particle-specs.dat', where the option is '-s' and the associated value is
    'particle-specs.dat'
    """

    # Defaults
    configuration = {'fn-trajectory': 'trajectory.dat',
                     'fn-particle-specs': 'particle-specs.dat',
                     'temperature': 298.15,
                     'E0': 0.0,
                     'direction-E0': 'x'}

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
        if item == '-d-E0':
            # Direction of electric external field.
            counter += 1
            configuration['direction-E0'] = argv[counter]
            entry = configuration['direction-E0']
            logging.info(f"Direction external field: {entry}")
        if item == '-T':
            # Temperature
            counter += 1
            configuration['temperature'] = float(argv[counter])
            entry = configuration['temperature']
            logging.info(f'Temperature: {entry}.')
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
        counter += 1

    return configuration
