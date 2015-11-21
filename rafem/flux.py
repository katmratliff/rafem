#! /usr/local/bin/python

import numpy as np

from avulsion_utils import get_link_lengths


def calc_qs(nu, riv_i, riv_j, n, dx, dy, dt):
    """Calculate sediment flux at river mouth."""
    ds = get_link_lengths((riv_i[-2:], riv_j[-2:]), dx=dx, dy=dy)
    dz = n[riv_i[-1], riv_j[-1]] - n[riv_i[-2], riv_j[-2]]

    return - nu * dz / ds
    #sed_flux = - (nu * dt) * ((n[riv_i[-1], riv_j[-1]] -
    #                           n[riv_i[-2], riv_j[-2]]) / dist)
