# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 20:19:37 2014

@author: kmratliff
"""
import numpy as np


def elev_change(imax, jmax, current_SL, n, riv_x, riv_y, ch_depth, dx, dy):

    rc_flag = np.zeros((imax, jmax))
    
    # creates river course flag array
    for i in range(len(riv_x)):
        
        rc_flag[riv_x[i]/dx][riv_y[i]/dy] = 1
        
        i += 1
    
    # changes elevation of last river course cell according to sea level change
    n[riv_x[-1]/dx][riv_y[-1]/dy] = current_SL - ch_depth

    # raises cell elevation to sea level if it is below
    for i in range(imax):
        for j in range(jmax):

            if (n[i][j] < current_SL and rc_flag[i][j] == 0):
                n[i][j] = current_SL

            j += 1
        i += 1
        
        # Need to somehow change above treatment of cells... they don't need to be
        # raised to sea level anymore. 
        
#    # raises elevation of whole inlet row
#    for I in range(jmax):
#        n[0][I] = n[0][I] + (IRR)
#        I = I + 1

    return n, rc_flag
