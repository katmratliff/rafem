# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 21:22:00 2015

@author: kmratliff
"""
import numpy as np


def cut_init(riv_i, riv_j, n, init_cut):
    n[riv_i, riv_j] -= init_cut

    return n


def cut_new(riv_i, riv_j, n, length_new, current_SL, a, ch_depth):
    """Set elevations of new portions of a profile."""

    # new last river cell = SL - channel depth
    n[riv_i[-1], riv_j[-1]] = current_SL - ch_depth
    
    if sum(length_new) > 0:
        n_new_segments = len(length_new)

        i0, j0 = riv_i[-n_new_segments], riv_j[-n_new_segments]

        # calculate slope of new stretch of river
        new_slope = (n[i0, j0] - (current_SL - ch_depth)) / sum(length_new)

        elevation_at_avulsion = n[riv_i[a], riv_j[a]]
        new_i, new_j = riv_i[-n_new_segments:], riv_j[-n_new_segments:]

        n[new_i, new_j] = elevation_at_avulsion - new_slope * length_new.cumsum()

    return n
