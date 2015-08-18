#! /usr/local/bin/python
import numpy as np

from avulsion_utils import is_diagonal_neighbor
from avulsion_utils import get_channel_distance


def solve_second_derivative(x, y):
    """
    Examples
    --------
    >>> import numpy as np
    >>> x = np.array([2., 3., 4.])
    >>> y = np.array([4., 9., 16.])
    >>> solve_second_derivative(x, y)
    """
    x2_minus_x1 = x[1:-1] - x[:-2]
    x3_minus_x2 = x[2:] - x[1:-1]
    x3_minus_x1 = x[2:] - x[:-2]

    return 2 * (y[:-2] / (x2_minus_x1 * x3_minus_x1) -
                y[1:-1] / (x3_minus_x2 * x2_minus_x1) +
                y[2:] / (x3_minus_x2 * x3_minus_x1))


def smooth_rc(dx, dy, nu, dt, riv_i, riv_j, n):
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
    n_river = n[riv_i, riv_j]
    s_river = get_channel_distance((riv_i, riv_j), dx=dx, dy=dy)

    dn_rc = (nu * dt) * solve_second_derivative(s_river, n_river)

    # NOTE: Divide by dx to match the old way, but I don't think this is
    # correct.
    dn_rc /= dx

    n[riv_i[1:-1], riv_j[1:-1]] += dn_rc

    return


# this function uses a linear diffusion equation (e.g. Paola 2000, Jerolmack
# and Paola 2007) to compute elevation change along the river course
def smooth_rc_old(dx, dy, nu, dt, riv_i, riv_j, n):
    # elevation change along river course due to diffusional smoothing
    for c in xrange(1, len(riv_i) - 1):
        n_prev = n[riv_i[c - 1], riv_j[c - 1]]
        n_cur = n[riv_i[c], riv_j[c]]
        n_next = n[riv_i[c + 1], riv_j[c + 1]]

        dwnst_dn = (n_next - n_cur) / dx
        upst_dn = (n_cur - n_prev) / dx

        if is_diagonal_neighbor((riv_i[c], riv_j[c]), (riv_i[c + 1], riv_j[c + 1])):
            dwnst_dn /= np.sqrt(2.)

        if is_diagonal_neighbor((riv_i[c], riv_j[c]), (riv_i[c - 1], riv_j[c - 1])):
            upst_dn /= np.sqrt(2.)

        dn_rc = (nu * dt) / (dx ** 2.) * (dwnst_dn - upst_dn)

        n[riv_i[c], riv_j[c]] += dn_rc

    return n
