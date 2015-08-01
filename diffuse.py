#! /usr/local/bin/python
import numpy as np


def is_diagonal_neighbor(sub0, sub1):
    return sub0[0] != sub1[0] and sub0[1] != sub1[1]


# this function uses a linear diffusion equation (e.g. Paola 2000, Jerolmack
# and Paola 2007) to compute elevation change along the river course
def smooth_rc(dx, dy, nu, dt, riv_i, riv_j, n):

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
