# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 21:22:00 2015

@author: kmratliff
"""
import numpy as np
import math

from avulsion_utils import get_link_lengths
from avulsion_utils import find_beach_length


def cut_init(riv_i, riv_j, n, init_cut):
    n[riv_i, riv_j] -= init_cut

    return n


def cut_new(riv_i, riv_j, n, sea_level, ch_depth, slope, dx=1., dy=1.):
    """Set elevations of new portions of a profile."""

    # downcut last river cell by a channel depth
    n[riv_i[-1], riv_j[-1]] -= ch_depth

    last_elev = n[riv_i[-1], riv_j[-1]] + ch_depth - sea_level
    max_cell_h = math.tan(slope) * dx

    beach_len = (last_elev / max_cell_h)

    ### NEED TO UPDATE INPUTS SO WE CAN USE BEACH LENGTH
    ### (JUST PASS PRIOR RIVER CELL, TOO) ###
    # beach_len = find_beach_length(n, (riv_i[-2], riv_j[-2]),
    #                               (riv_i[-1], riv_j[-1]), sea_level,
    #                               ch_depth, slope, dx=dx, dy=dy)
    
    if riv_i.size > 1:
        lengths = get_link_lengths((riv_i, riv_j), dx=dx, dy=dy)
        lengths[-1] = np.divide(lengths[-1], 2) + beach_len

        i0, j0 = riv_i[1], riv_j[1]
        z0 = n[riv_i[0], riv_j[0]]

        # calculate slope of new stretch of river
        slope = (n[i0, j0] - n[riv_i[-1], riv_j[-1]]) / lengths.sum()

        n[riv_i[1:], riv_j[1:]] = z0 - slope * lengths.cumsum()
