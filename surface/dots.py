"""Dotted surface"""

import numpy as np
import logging


def read(fn: str) -> np.ndarray:
    """Reads points of a dotted surface representation.
    :param fn Input file name.
    :return number of dots (n_dots), Array of points with shape = (n_dots, 3)
    """
    f = open(fn, 'r');

    # First line.
    line = f.readline()
    items = line.split()
    n_dots = int(items[0])

    # Read the dots
    dots = np.ndarray(shape=(n_dots, 3))
    counter = 0
    for _ in np.arange(0, n_dots):
        line = f.readline()
        items = line.split()
        dot = np.array([float(items[0]), float(items[1]), float(items[2])])
        dots[counter] = dot
        counter += 1

    f.close()
    logging.info(f'Read {n_dots} surface dots.')

    return n_dots, dots
