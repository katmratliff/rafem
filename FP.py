# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 08:56:09 2014

@author: kmratliff
"""
import numpy as np


def dep_blanket(current_SL, blanket_rate, n, riv_i, riv_j, ch_depth):
    depo_flag = np.ones(n.shape, dtype=np.int)

    # don't deposit in the river course
    depo_flag[riv_i, riv_j] = 0

    # don't deposit if cell <= sea level
    # ASK BRAD ABOUT THIS (no cell is below SL right now)
    # but it might be when couple?? need to figure out if keeping up
    # with elevation RELATIVE to sea level or actual elevation??
    depo_flag[n <= current_SL] = 0

    # if cell elevation is above bankfull elev, don't deposit
    bankfull_elevation = n[riv_i, riv_j] + ch_depth
    for row in riv_i:
        depo_flag[row, n[row] > bankfull_elevation[row]] = 0

    # don't deposit on first two rows b/c inlet rise rate does that
    depo_flag[:2, :] = 0

    # deposit "blanket" deposition on qualified cells
    n[depo_flag == 1] += blanket_rate

    dn_fp = depo_flag * blanket_rate

    return n, dn_fp


def distance_to_river(y, y0):
    return np.absolute(y - y0)


def within_wetland(y, riv_ind, wetland_width=0.):
    dy = distance_to_river(y, y[riv_ind])
    is_wetland = dy < wetland_width
    is_wetland[riv_ind] = False
    return is_wetland


def wetlands(current_SL, WL_Z, wetland_width, n, riv_i, riv_j, x, y, dn_fp):
    depo_wetland = np.zeros(n.shape, dtype=np.int)

    for row, col in zip(riv_i, riv_j):
        dist = within_wetland(y[row], col, wetland_width=wetland_width)
        elev = n[row] < current_SL + WL_Z

        cols = dist & elev & depo_wetland

        before = n[row, cols].copy()
        n[row, cols] = current_SL + WL_Z
        wetland_dep = n[row, cols] - before

        dn_fp[row, cols] += wetland_dep
        depo_wetland[row, cols] == 1

    return n, dn_fp


def dep_splay(riv_i, riv_j, new_riv_i, new_riv_j,
              ch_depth, n, a, dn_fp, splay_type, splay_dep):
    """
    USE DEPTH-DEPENDENCE IN THE FUTURE
        SE1 = (n[riv_x[a]/dx][riv_y[a]/dy] + ch_depth - \
                n[new_riv_x[a]/dx][new_riv_y[a]/dy]) / ch_depth

    ADD LARGER SPLAY TYPE IN FUTURE?
    (could be first two failed river cells + surrounding)

    This could possibly be improved by comparing to find nearest beach routine (CEM)
    or using some sort of search radius 
    """
    i_new, j_new = new_riv_i[a], new_riv_j[a]

    if splay_type == 1:  # splay deposition just at first failed river cell

        # deposit at failed avulsion river cell
        if i_new != riv_i[a] or j_new != riv_j[a]:
            n[i_new, j_new] += splay_dep

        # record deposition in dn_fp
        dn_fp[i_new, j_new] += splay_dep

    if splay_type == 2:     # splay deposition at first failed river cell
                            # and the adjacent cells

        depo_flag2 = np.zeros(n.shape, dtype=np.int)

        if i_new == n.shape[0]:
            depo_flag2[i_new, j_new - 1] = 1    # left side
            depo_flag2[i_new, j_new + 1] = 1    # right side
            depo_flag2[i_new - 1, j_new - 1] = 1  # Left u.s.
            depo_flag2[i_new - 1, j_new] = 1  # center u.s.
            depo_flag2[i_new - 1, j_new + 1] = 1  # right u.s.
        elif j_new == 0:
            depo_flag2[i_new, j_new] = 1  # failed river cell
            depo_flag2[i_new + 1, j_new] = 1  # center d.s.
            depo_flag2[i_new + 1, j_new + 1] = 1  # right d.s.
            depo_flag2[i_new, j_new + 1] = 1    # right side
            depo_flag2[i_new - 1, j_new] = 1  # center u.s.
            depo_flag2[i_new - 1, j_new + 1] = 1  # right u.s.
        elif j_new == n.shape[1]:
            depo_flag2[i_new, j_new] = 1  # failed river cell
            depo_flag2[j_new + 1, j_new - 1] = 1  # Left d.s.
            depo_flag2[i_new + 1, j_new] = 1  # center d.s.
            depo_flag2[i_new, j_new - 1] = 1    # left side
            depo_flag2[i_new - 1, j_new - 1] = 1  # Left u.s.
            depo_flag2[i_new - 1, j_new] = 1  # center u.s.
        else:         
            depo_flag2[i_new, j_new] = 1  # failed river cell
            depo_flag2[i_new + 1, j_new - 1] = 1  # Left d.s.
            depo_flag2[i_new + 1, j_new] = 1  # center d.s.
            depo_flag2[i_new + 1, j_new + 1] = 1  # right d.s.
            depo_flag2[i_new, j_new - 1] = 1    # left side
            depo_flag2[i_new, j_new + 1] = 1    # right side
            depo_flag2[i_new - 1, j_new - 1] = 1  # Left u.s.
            depo_flag2[i_new - 1, j_new] = 1  # center u.s.
            depo_flag2[i_new - 1, j_new + 1] = 1  # right u.s.

        # no splay deposition in river channel
        depo_flag2[riv_i, riv_j] = 0

            i += 1

        # deposit splay sediment on flagged cells
        n[depo_flag2 == 1] += splay_dep
        dn_fp[depo_flag2 == 1] += splay_dep

                j += 1
            i += 1

    return n, dn_fp
