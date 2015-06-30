# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 08:56:09 2014

@author: kmratliff
"""
import numpy as np


def dep_blanket(dy, dx, imax, jmax, current_SL, blanket_rate, n,
                riv_x, riv_y, ch_depth):

    depo_flag = np.ones((imax, jmax))

    # don't deposit in the river course
    for k in range(len(riv_x)):

        depo_flag[riv_x[k]/dx][riv_y[k]/dy] = 0

        k += 1

    for i in range(imax):
        for j in range(jmax):

            # don't deposit if cell <= sea level
            """
            ASK BRAD ABOUT THIS (no cell is below SL right now)
            but it might be when couple?? need to figure out if keeping up
            with elevation RELATIVE to sea level or actual elevation??
            """
            if n[i][j] <= current_SL:

                depo_flag[i][j] = 0
            
            j += 1
        i += 1

    # if cell elevation is above bankfull elev, don't deposit
    for i in range(len(riv_x)):
        for j in range(jmax):
            
            if (n[riv_x[i]/dx][j] >= (n[riv_x[i]/dx][riv_y[i]/dy]
                                      + ch_depth)):

                depo_flag[riv_x[i]/dx][j] = 0

            j += 1
        i += 1
        
    # don't deposit on first two rows b/c inlet rise rate does that
    depo_flag[0][:] = 0
    depo_flag[1][:] = 0

    for i in range(imax):
        for j in range(jmax):
            # deposit "blanket" deposition on qualified cells
            if depo_flag[i][j] == 1:

                n[i][j] = n[i][j] + blanket_rate

            j += 1
        i += 1

    dn_fp = depo_flag * blanket_rate

    return n, dn_fp

def wetlands(dx, dy, imax, jmax, current_SL, WL_Z, WL_dist, n, riv_x, riv_y,
             x, y, dn_fp):
    
    depo_wetland = np.zeros((imax,jmax))

    for i in range(len(riv_x)):
        for j in range(jmax):
            
            # determine cells within specified wetland distance from river
            if (np.absolute(y[riv_x[i]/dx][riv_y[i]/dy] - y[riv_x[i]/dx][j]) \
                <= (WL_dist * dy) and (y[riv_x[i]/dx][riv_y[i]/dy] - \
                y[riv_x[i]/dx][j]) != 0):
                dist = 1
            else:
                dist = 0

            # determine if a cell is below wetland elevation + current SL
            if n[riv_x[i]/dx][j] < (current_SL + WL_Z):
                elev = 1
            else:
                elev = 0
            
            if (dist == 1 and elev == 1 and depo_wetland[riv_x[i]/dx][j] == 0):
                
                before = n[riv_x[i]/dx][j]          
                n[riv_x[i]/dx][j] = current_SL + WL_Z
                
                # record "wetland deposition"
                wetland_dep = n[riv_x[i]/dx][j] - before
                dn_fp[riv_x[i]/dx][j] = dn_fp[riv_x[i]/dx][j] + wetland_dep
            
                # depo flag so deposition at cell isn't recorded more than once
                depo_wetland[riv_x[i]/dx][j] == 1

            j += 1
        i += 1
    
    return n, dn_fp

def dep_splay(dy, dx, imax, jmax, riv_x, riv_y, new_riv_x, new_riv_y,
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
    if splay_type == 1:  # splay deposition just at first failed river cell

        # deposit at failed avulsion river cell
        if new_riv_x[a] != riv_x[a] or new_riv_y[a] != riv_y[a]:
            n[new_riv_x[a]/dx][new_riv_y[a]/dy] = \
                n[new_riv_x[a]/dx][new_riv_y[a]/dy] + splay_dep

        # record deposition in dn_fp
        dn_fp[new_riv_x[a]/dx][new_riv_y[a]/dy] = \
            dn_fp[new_riv_x[a]/dx][new_riv_y[a]/dy] + splay_dep

    if splay_type == 2:     # splay deposition at first failed river cell
                            # and the adjacent cells

        depo_flag2 = np.zeros((imax, jmax))

        if new_riv_x[a]/dx == imax:
            
            depo_flag2[new_riv_x[a]/dx][(new_riv_y[a]/dy)-1] = 1    # left side
            depo_flag2[new_riv_x[a]/dx][(new_riv_y[a]/dy)+1] = 1    # right side
            depo_flag2[(new_riv_x[a]/dx)-1][(new_riv_y[a]/dy)-1] = 1  # Left u.s.
            depo_flag2[(new_riv_x[a]/dx)-1][new_riv_y[a]/dy] = 1  # center u.s.
            depo_flag2[(new_riv_x[a]/dx)-1][(new_riv_y[a]/dy)+1] = 1  # right u.s.
            
        elif new_riv_y[a]/dy == 0:
            
            depo_flag2[new_riv_x[a]/dx][new_riv_y[a]/dy] = 1  # failed river cell
            depo_flag2[(new_riv_x[a]/dx)+1][new_riv_y[a]/dy] = 1  # center d.s.
            depo_flag2[(new_riv_x[a]/dx)+1][(new_riv_y[a]/dy)+1] = 1  # right d.s.
            depo_flag2[new_riv_x[a]/dx][(new_riv_y[a]/dy)+1] = 1    # right side
            depo_flag2[(new_riv_x[a]/dx)-1][new_riv_y[a]/dy] = 1  # center u.s.
            depo_flag2[(new_riv_x[a]/dx)-1][(new_riv_y[a]/dy)+1] = 1  # right u.s.
            
        elif new_riv_y[a]/dy == jmax:
            
            depo_flag2[new_riv_x[a]/dx][new_riv_y[a]/dy] = 1  # failed river cell
            depo_flag2[(new_riv_x[a]/dx)+1][(new_riv_y[a]/dy)-1] = 1  # Left d.s.
            depo_flag2[(new_riv_x[a]/dx)+1][new_riv_y[a]/dy] = 1  # center d.s.
            depo_flag2[new_riv_x[a]/dx][(new_riv_y[a]/dy)-1] = 1    # left side
            depo_flag2[(new_riv_x[a]/dx)-1][(new_riv_y[a]/dy)-1] = 1  # Left u.s.
            depo_flag2[(new_riv_x[a]/dx)-1][new_riv_y[a]/dy] = 1  # center u.s.
        
        else:         
        
            depo_flag2[new_riv_x[a]/dx][new_riv_y[a]/dy] = 1  # failed river cell
            depo_flag2[(new_riv_x[a]/dx)+1][(new_riv_y[a]/dy)-1] = 1  # Left d.s.
            depo_flag2[(new_riv_x[a]/dx)+1][new_riv_y[a]/dy] = 1  # center d.s.
            depo_flag2[(new_riv_x[a]/dx)+1][(new_riv_y[a]/dy)+1] = 1  # right d.s.
            depo_flag2[new_riv_x[a]/dx][(new_riv_y[a]/dy)-1] = 1    # left side
            depo_flag2[new_riv_x[a]/dx][(new_riv_y[a]/dy)+1] = 1    # right side
            depo_flag2[(new_riv_x[a]/dx)-1][(new_riv_y[a]/dy)-1] = 1  # Left u.s.
            depo_flag2[(new_riv_x[a]/dx)-1][new_riv_y[a]/dy] = 1  # center u.s.
            depo_flag2[(new_riv_x[a]/dx)-1][(new_riv_y[a]/dy)+1] = 1  # right u.s.

        # no splay deposition in river channel
        for i in range(len(riv_x)):

            depo_flag2[riv_x[i]/dx][riv_y[i]/dy] = 0

            i += 1

        # deposit splay sediment on flagged cells
        for i in range(imax):
            for j in range(jmax):
                if depo_flag2[i][j] == 1:

                    n[i][j] = n[i][j] + splay_dep
                    dn_fp[i][j] = dn_fp[i][j] + splay_dep

                j += 1
            i += 1

    return n, dn_fp
