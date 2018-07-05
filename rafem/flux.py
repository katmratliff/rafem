#! /usr/local/bin/python

import numpy as np

from .avulsion_utils import get_link_lengths, find_beach_length_riv_cell


def calc_qs(nu, riv_i, riv_j, n, sea_level, ch_depth, dx, dy, dt, slope):
    """Calculate sediment flux at river mouth."""

    beach_len = find_beach_length_riv_cell(n, (riv_i[-2], riv_j[-2]),
                                  (riv_i[-1], riv_j[-1]), sea_level,
                                  ch_depth, slope, dx=dx, dy=dy)

    ds = get_link_lengths((riv_i[-2:], riv_j[-2:]), dx=dx, dy=dy)
    ds[-1] += beach_len
    dz = (sea_level - ch_depth) - n[riv_i[-2], riv_j[-2]]

    return - nu * dz / ds

    # ds = get_link_lengths((riv_i[-3:-1], riv_j[-3:-1]), dx=dx, dy=dy)
    
    # dz = n[riv_i[-2], riv_j[-2]] - n[riv_i[-3], riv_j[-3]]

    # return - nu * dz / ds




