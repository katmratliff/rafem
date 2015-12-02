#! /usr/local/bin/python
import numpy as np

from avulsion_utils import is_diagonal_neighbor
from avulsion_utils import get_channel_distance


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


def smooth_rc(dx, dy, nu, dt, riv_i, riv_j, n, sea_level):
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
    # NOTE: Divide by dx to match the old way, but I don't think this is
    # correct.
    # nu /= dx
    # KMR 8/24/15: don't need to divide by dx anymore, diffusion coeff
    # should be fixed with new calculation

    beach_len = n[riv_i[-1], riv_j[-1]] - sea_level

    if beach_len > 1:
        n_river = n[riv_i, riv_j]
        s_river = get_channel_distance((riv_i, riv_j), dx=dx, dy=dy)

        dn_rc = (nu * dt) * solve_second_derivative(s_river, n_river)

        n[riv_i[1:-1], riv_j[1:-1]] += dn_rc

    else:
        n_river = n[riv_i[:-1], riv_j[:-1]]
        s_river = get_channel_distance((riv_i[:-1], riv_j[:-1]), dx=dx, dy=dy)

        dn_rc = (nu * dt) * solve_second_derivative(s_river, n_river)

        n[riv_i[1:-2], riv_j[1:-2]] += dn_rc

        n_river_last = n[riv_i[-3:], riv_j[-3:]]
        s_river_last = get_channel_distance((riv_i[-3:], riv_j[-3:]), dx=dx, dy=dy)
        np.divide(s_river_last[-1], 2)
        s_river_last[-1] += beach_len

        dn_rc_last = (nu * dt) * solve_second_derivative(s_river_last,
                                                         n_river_last)

        n[riv_i[-2], riv_j[-2]] += dn_rc_last

    return


# this function uses a linear diffusion equation (e.g. Paola 2000, Jerolmack
# and Paola 2007) to compute elevation change along the river course
def smooth_rc_old(dx, dy, nu, dt, riv_i, riv_j, n):
    # elevation change along river course due to diffusional smoothing
    for c in xrange(1, len(riv_i) - 1):
        n_prev = n[riv_i[c - 1], riv_j[c - 1]]
        n_cur = n[riv_i[c], riv_j[c]]
        n_next = n[riv_i[c + 1], riv_j[c + 1]]

        dwnst_dx, upst_dx = dx, dx
        if is_diagonal_neighbor((riv_i[c], riv_j[c]), (riv_i[c + 1], riv_j[c + 1])):
            dwnst_dx *= np.sqrt(2.)

        if is_diagonal_neighbor((riv_i[c], riv_j[c]), (riv_i[c - 1], riv_j[c - 1])):
            upst_dx *= np.sqrt(2.)

        dwnst_dn = (n_next - n_cur) / dwnst_dx
        upst_dn = (n_cur - n_prev) / upst_dx
        mean_dx = (dwnst_dx + upst_dx) * .5

        # NOTE: This is the old way but, I think, is incorrect. For
        # non-uniform spacing of points you need to divide by the mean spacing.
        #dn_rc = (nu * dt) / (dx ** 2.) * (dwnst_dn - upst_dn)

        # This properly solves the second derivative with unequal spacing in x.
        dn_rc = (nu / dx * dt) * (dwnst_dn - upst_dn) / mean_dx

        n[riv_i[c], riv_j[c]] += dn_rc

    return n
