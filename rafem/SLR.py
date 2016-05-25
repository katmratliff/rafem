#! /usr/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np


def elev_change(current_SL, n, riv_i, riv_j, ch_depth, SLRR):
    """Raise elevations to sea level if below sea level.

    Set elevations of cells below sea level to sea level unless they are
    within the river channel.

    Parameters
    ----------
    current_SL : float
        Sea level.
    n : ndarray
        Array of elevations.
    riv_i : ndarray of int
        Row indices into *n* for river cells.
    riv_j : ndarray of int
        Column indices into *n* for river cells.
    ch_depth : float
        Channel depth.
    """
    # changes elevation of last river course cell according to sea level change
    #n[riv_i[-1], riv_j[-1]] = current_SL - ch_depth

    n[riv_i[-1], riv_j[-1]] += SLRR

    channel_elevations = n[riv_i, riv_j].copy()

    # raises cell elevation to sea level if it is below
    # n[n < current_SL] = current_SL
    #channel_elevations[channel_elevations < current_SL] = current_SL

    n[riv_i, riv_j] = channel_elevations

    # Need to somehow change above treatment of cells... they don't need to be
    # raised to sea level anymore. 
        
#    # raises elevation of whole inlet row
#    for I in range(jmax):
#        n[0][I] = n[0][I] + (IRR)
#        I = I + 1
