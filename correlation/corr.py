import numpy as np
import tidynamics as td


def acf(data, n_use=50):
    """Computes the auto correlation function (acf) of a time dependent quantity or vector

    Parameters
    ----------
        data : ndarray
            An (N,M) numpy array of dimension 2. The first column should
            represent time, while the remaining ones are quantity or
            vector components values.
        n_use : int
            Only the first n_use elements of the ACF are returned. Default is 50.

    Returns
    -------
        acf(data) : tuple
            An tuple of length 3 with time, zeros, and acf, where acf is
            normalized, each of length n_use.
    """
    time = data[:, 0]  # Each row, first column.
    if time.size < n_use:
        msg = 'Cannot select first ' + str(n_use) + ' elements from ACF of length ' + str(time.size) + '.'
        raise Exception(msg)
    time0 = time[0]
    time -= time0
    v = data[:, 1:]    # Each row, columns 1 and up.
    acf = td.acf(v)
    acf0 = acf[0]
    nacf = acf / acf0
    pacf = nacf[0:n_use]
    ptime = time[0:n_use]
    zeros = np.zeros(n_use)
    return ptime, zeros, pacf
