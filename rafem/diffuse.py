#! /usr/local/bin/python

import numpy as np

from .avulsion_utils import (
    is_diagonal_neighbor,
    get_channel_distance,
    find_beach_length_riv_cell,
)


def solve_second_derivative(x, y):
    """Solve the second derivative of y with respect to x.

    Use finite difference to solve the second derivative of *y* with respect
    to *x* where *x* can be unevenly spaced.

    Examples
    --------
    Values can be evenly spaced.

    >>> import numpy as np
    >>> x = np.array([2., 3., 4.])
    >>> y = np.array([4., 9., 16.])
    >>> solve_second_derivative(x, y)
    array([ 2.])

    Values are unevenly spaced.

    >>> x = np.array([2., 3., 5.])
    >>> y = np.array([4., 9., 25.])
    >>> solve_second_derivative(x, y)
    array([ 2.])
    """
    x2_minus_x1 = x[1:-1] - x[:-2]
    x3_minus_x2 = x[2:] - x[1:-1]
    x3_minus_x1 = x[2:] - x[:-2]

    return 2 * (y[:-2] / (x2_minus_x1 * x3_minus_x1) -
                y[1:-1] / (x3_minus_x2 * x2_minus_x1) +
                y[2:] / (x3_minus_x2 * x3_minus_x1))


def smooth_rc(dx, dy, nu, dt, ch_depth, riv_i, riv_j, n, sea_level, slope):
    """Smooth river channel elevations using the diffusion equation.

    Parameters
    ----------
    dx : float
        Spacing of grid columns.
    dy : float
        Spacing of grid rows.
    nu : float
        Diffusion coefficient.
    dt : float
        Time step (in seconds).
    riv_i : ndarray
        Row indices for the river path.
    riv_j : ndarray
        Column indices for the river path.
    n : ndarray
        2D array of grid elevations.
    """

    beach_len = find_beach_length_riv_cell(n, (riv_i[-2], riv_j[-2]),
                                  (riv_i[-1], riv_j[-1]), sea_level,
                                  ch_depth, slope, dx=dx, dy=dy)


    n_river = n[riv_i, riv_j]
    n_river[-1] = sea_level - ch_depth
    s_river = get_channel_distance((riv_i, riv_j), dx=dx, dy=dy)
    s_river[-1] += beach_len

    dn_rc = (nu * dt) * solve_second_derivative(s_river, n_river)

    n[riv_i[1:-1], riv_j[1:-1]] += dn_rc

    return dn_rc

def calc_crevasse_dep(dx, dy, nu, dt, ch_depth, riv_i, riv_j, n,
                      sea_level, slope, loc):
    """Calculate crevasse splay deposition rate."""

    beach_len = find_beach_length_riv_cell(n, (riv_i[-2], riv_j[-2]),
                                  (riv_i[-1], riv_j[-1]), sea_level,
                                  ch_depth, slope, dx=dx, dy=dy)

    n_river = n[riv_i, riv_j]
    n_river[-1] = sea_level - ch_depth
    s_river = get_channel_distance((riv_i, riv_j), dx=dx, dy=dy)
    s_river[-1] += beach_len

    dn_rc = (nu * dt) * solve_second_derivative(s_river, n_river)

    # average deposition in 3 cells above 'breach' to find
    # crevasse splay deposition rate
    if len(dn_rc[:loc]) >= 3:
        splay_dep = np.average(dn_rc[loc-3:loc])
    else: splay_dep = dn_rc[loc-1]

    if splay_dep < 0:
        splay_dep = 0

    return splay_dep

