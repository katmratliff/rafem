#! /usr/local/bin/python

import numpy as np
import pudb


def linear_subsidence(n, riv_i, riv_j, ch_depth, sub_rate, sub_start, SL):
    """ Subside cells in rows beyond start location by linear rate.
    Updated version could be to subside at an increasing rate towards ocean. """

    pu.db

    subside_cells = np.zeros_like(n)
    subside_cells[sub_start:,:] = 1

    subaerial_elev = n
    subaerial_elev[riv_i, riv_j] += ch_depth

    subside_cells[n <= SL] = 0

    n[subside_cells == 1] -= sub_rate