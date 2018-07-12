#! /usr/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import math

from .avulsion_utils import get_link_lengths, find_new_beach_length


def cut_init(riv_i, riv_j, n, init_cut):
    n[riv_i, riv_j] -= init_cut

    return n


def cut_new(riv_i, riv_j, n, sea_level, ch_depth, dx=1., dy=1.):
    """Set elevations of new portions of a profile."""
    
    if riv_i.size > 1:
        beach_len = find_new_beach_length(n, (riv_i[-2], riv_j[-2]),
                                      (riv_i[-1], riv_j[-1]), sea_level,
                                      dx=dx, dy=dy)

        lengths = get_link_lengths((riv_i, riv_j), dx=dx, dy=dy)
        lengths[-1] += beach_len 

        i0, j0 = riv_i[0], riv_j[0]
        z0 = n[riv_i[0], riv_j[0]]

        # calculate slope of new stretch of river
        # slope = (n[i0, j0] - n[riv_i[-1], riv_j[-1]]) / lengths.sum()
        slope = ((n[i0, j0] - (n[riv_i[-1], riv_j[-1]] - ch_depth))
                / lengths.sum())

        n[riv_i[1:], riv_j[1:]] = z0 - slope * lengths.cumsum()


def cut_local(riv_i, riv_j, n, dx=1., dy=1.):
    """Set elevations of new portions of a profile after local avulsion."""
    
    if riv_i.size > 1:

        lengths = get_link_lengths((riv_i, riv_j), dx=dx, dy=dy) 

        i0, j0 = riv_i[0], riv_j[0]
        z0 = n[riv_i[0], riv_j[0]]

        # calculate slope of new stretch of river
        slope = (n[i0, j0] - n[riv_i[-1], riv_j[-1]]) / lengths.sum()

        n[riv_i[1:], riv_j[1:]] = z0 - slope * lengths.cumsum()

