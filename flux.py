#! /usr/local/bin/python

import numpy as np


def is_diagonal_neighbor(sub0, sub1):
    return sub0[0] != sub1[0] and sub0[1] != sub1[1]


def calc_qs(nu, riv_i, riv_j, n, dx, dy, dt):

    sed_flux = 0
    dist = 0

    if is_diagonal_neighbor((riv_i[-1], riv_j[-1]), (riv_i[-2], riv_j[-2])):
        dist = np.sqrt(2.)
    else:
        dist = 1.

    sed_flux = - (nu * dt) * ((n[riv_i[-1], riv_j[-1]] -
                               n[riv_i[-2], riv_j[-2]]) / dist)

    return sed_flux
