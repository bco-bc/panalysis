"""Triangulated surface
"""

import numpy as np
import logging

def read(fn: str) -> [int, np.ndarray]:
    """Reads a triangulated surface. An edge consists of a start and an end point (vertex).
    :param fn Input file name.
    :return Number of edges (n_edges), Array of edges with shape = (n_edges, 6) where the first (second) set
    of 3 values is the begin (end) point of the edge.
    """
    f = open(fn, 'r')

    # First line
    line = f.readline()
    items = line.split()


    f.close()
    # logging.info(f'Read {n_edges} edges.')
    return 0
