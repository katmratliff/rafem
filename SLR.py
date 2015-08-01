# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 20:19:37 2014

@author: kmratliff
"""
import numpy as np


def elev_change(current_SL, n, riv_i, riv_j, ch_depth):
    # creates river course flag array
    rc_flag = np.zeros(n.shape, dtype=np.int)
    rc_flag[riv_i, riv_j] = 1
    
    # changes elevation of last river course cell according to sea level change
    n[riv_i[-1], riv_j[-1]] = current_SL - ch_depth

    # raises cell elevation to sea level if it is below
    n[(n < current_SL) & (rc_flag == 0)] = current_SL

    # Need to somehow change above treatment of cells... they don't need to be
    # raised to sea level anymore. 
        
#    # raises elevation of whole inlet row
#    for I in range(jmax):
#        n[0][I] = n[0][I] + (IRR)
#        I = I + 1

    return n
