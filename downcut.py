# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 21:22:00 2015

@author: kmratliff
"""

def cut_init(dx, dy, riv_x, riv_y, n, init_cut):
    
#    n[riv_x[:]/dx][riv_y[:]/dy] = n[riv_x[:]/dx][riv_y[:]/dy] - init_cut
    
    for i in range(len(riv_x)):
        
        n[riv_x[i]/dx][riv_y[i]/dy] = n[riv_x[i]/dx][riv_y[i]/dy] - init_cut
        
#        if n[riv_x[i]/dx][riv_y[i]/dy] < Initial_SL:
#            n[riv_x[i]/dx][riv_y[i]/dy] = Initial_SL
        
        i += 1
            
    return n

def cut_new(dx, dy, riv_x, riv_y, n, length_new, current_SL, a, ch_depth):

    # new last river cell = SL - channel depth
    n[riv_x[-1]/dx][riv_y[-1]/dy] = current_SL - ch_depth
    
    if sum(length_new) > 0:

        # calculate slope of new stretch of river
        new_slope = ((n[riv_x[-len(length_new)]/dx][riv_y[-len(length_new)]/dy]
                        - (current_SL - ch_depth)) / sum(length_new))
        
        for k in range(1,len(riv_x[a:])-1):

            n[riv_x[a+k]/dx][riv_y[a+k]/dy] = (n[riv_x[a]/dx][riv_y[a]/dy] - 
                                                new_slope * sum(length_new[:k]))
            
            k += 1
            
    return n