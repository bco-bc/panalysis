"""Edges polyhedron
"""

import numpy as np
import logging


def read(fn: str) -> [int, np.ndarray]:
    """Reads edges of a polyhedron
    :param fn Input fine name
    :return Number of edges (n_edges), Array of edges with shape = (n_edges, 6), with a start point and
    end point
    """
    f = open(fn, 'r')

    # First line
    line = f.readline()
    items = line.split()
    n_edges = int(items[0])

    # Read edges.
    edges = np.ndarray(shape=(n_edges, 6))
    counter = 0
    for _ in np.arange(0, n_edges):
        line = f.readline()
        items = line.split()
        edge = np.array([float(items[0]), float(items[1]), float(items[2]),
                         float(items[3]), float(items[4]), float(items[5])])
        edges[counter] = edge
        counter += 1

    f.close()
    logging.info(f'Read {n_edges} edges.')
    return n_edges, edges
