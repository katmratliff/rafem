#! /usr/local/bin/python
import numpy as np


def lowest_neighbor(n, i, j):
    """Find lowest neighbor value around a point."""
    if j == n.shape[1] - 1:
        di, dj  = np.array([0, 1, 1]), np.array([-1, -1, 0])
    elif j == 0:
        di, dj  = np.array([1, 1, 0]), np.array([0, 1, 1])
    else:
        di, dj = np.array([0, 1, 1, 1, 0]),  np.array([-1, -1, 0, 1, 1])

    lowest = np.argmin(n[i + di, j + dj])
    return i + di[lowest], j + dj[lowest]


def find_course(n, riv_i, riv_j, sea_level=None):
    # function to find the steepest descent route
    # note: this needs to be improved to remove potential bias that may occur
    # if two or more cells have the steepest descent elevation

    riv_len = 1
    while 1:
        if riv_i[riv_len - 1] == n.shape[0] - 1 or riv_len == riv_i.size:
            break
        if sea_level is not None and n[riv_i[riv_len - 1], riv_j[riv_len - 1]] <= sea_level:
            break
        riv_i[riv_len], riv_j[riv_len] = lowest_neighbor(n, riv_i[riv_len - 1], riv_j[riv_len - 1])
        riv_len += 1

    return riv_len
