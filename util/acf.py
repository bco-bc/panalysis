import tidynamics as td


def calculate(data, n_use=50):
    """Calculates the normalized auto util function (ACF) of a time dependent quantity or vector

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
        tuple
            An tuple with time and ACF, where ACF is
            normalized, each of length n_use.

    Raises
    ------
        ValueError
            If the value of n_use is larger than the length of ACF.
    """
    time = data[:, 0]  # Each row, first column.
    if time.size < n_use:
        msg = 'Cannot select first ' + str(n_use) + ' elements from ACF of length ' + str(time.size) + '.'
        raise ValueError(msg)
    time0 = time[0]
    time -= time0
    v = data[:, 1:]    # Each row, columns 1 and up.
    acf = td.acf(v)
    acf0 = acf[0]
    nacf = acf / acf0
    pacf = nacf[0:n_use]
    ptime = time[0:n_use]
    return ptime, pacf
