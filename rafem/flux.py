#! /usr/local/bin/python

import numpy as np

from avulsion_utils import get_link_lengths
from avulsion_utils import is_diagonal_neighbor


def calc_qs(nu, riv_i, riv_j, n, sea_level, ch_depth, dx, dy, dt):
    """Calculate sediment flux at river mouth."""

    beach_len = n[riv_i[-1], riv_j[-1]] + ch_depth - sea_level

    ds = get_link_lengths((riv_i[-2:], riv_j[-2:]), dx=dx, dy=dy)
    dz = n[riv_i[-1], riv_j[-1]] - n[riv_i[-2], riv_j[-2]]

    if beach_len >= 1:
    	flux = - nu * dz / ds

    else:
    	flux = - nu * dz / (ds/2 + beach_len)


    return - nu * dz / ds
    #sed_flux = - (nu * dt) * ((n[riv_i[-1], riv_j[-1]] -
    #                           n[riv_i[-2], riv_j[-2]]) / dist)
