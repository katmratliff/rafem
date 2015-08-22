#! /usr/local/bin/python

import steep_desc
import downcut
import FP
import numpy as np
import math


# determines if there is an avulsion along river course
def find_avulsion(dx, dy, imax, jmax, riv_x, riv_y, n, super_ratio, current_SL,
                  ch_depth, short_path, dn_fp, splay_type, splay_dep):

    loc = []
    SEL = np.zeros(len(riv_x))
    SER = np.zeros(len(riv_x))
    avulsion_type = 0
    length_new_sum = 0    
    length_old = 0

    for a in range(1, len(riv_x)):

        ch_Z = n[riv_x[a]/dx][riv_y[a]/dy] + ch_depth   # bankfull elev.
        LHS = n[riv_x[a]/dx][(riv_y[a]/dy)-1]
        RHS = n[riv_x[a]/dx][(riv_y[a]/dy)+1]

        # normalized superelevation ratio on left side
        SEL[a] = ((ch_Z - LHS) / ch_depth)

        # normalized superelevation ratio on right side
        SER[a] = ((ch_Z - RHS) / ch_depth)

        if SEL[a] >= super_ratio or SER[a] >= super_ratio:

            # if superelevation greater than trigger ratio, determine
            # length of new steepest descent path
            new_riv_x = riv_x[:a-1]
            new_riv_y = riv_y[:a-1]

            new_riv_x, new_riv_y = steep_desc.find_new_course(
                dx, dy, imax, jmax, n, new_riv_x, new_riv_y, current_SL)

            # if using the shortest path as an avulsion criterion, then
            # the lengths of the previous and newly calculated paths will
            # be compared
            if short_path == 1:
                
                # duplicates arrays so that length can be compared below
                test_new_x = new_riv_x[a:]
                test_new_y = new_riv_y[a:]
                test_old_x = riv_x[a:]
                test_old_y = riv_y[a:]
                length_new = []
                
                for c in range(len(test_new_x)-1):
                    
                    if (((test_new_x[c+1]/dx) - (test_new_x[c]/dx) == 0) and
                        (test_new_y[c+1]/dy) - (test_new_y[c]/dy) == -1):
                            
                            length_new.append(1)
                    
                    elif (((test_new_x[c+1]/dx) - (test_new_x[c]/dx) == 0)
                        and (test_new_y[c+1]/dy) - (test_new_y[c]/dy) == 1):
                    
                            length_new.append(1)
                            
                    elif (((test_new_x[c+1]/dx) - (test_new_x[c]/dx) == 1)
                        and (test_new_y[c+1]/dy) - (test_new_y[c]/dy) == 0):

                            length_new.append(1)
                    
                    elif (((test_new_x[c+1]/dx) - (test_new_x[c]/dx) == 1)
                        and (test_new_y[c+1]/dy) - (test_new_y[c]/dy) == -1):
                            
                            length_new.append(math.sqrt(2))

                    elif (((test_new_x[c+1]/dx) - (test_new_x[c]/dx) == 1)
                        and (test_new_y[c+1]/dy) - (test_new_y[c]/dy) == 1):
                            
                            length_new.append(math.sqrt(2))
                
                for b in range(len(test_old_x)-1):
                    
                    if (((test_old_x[b+1]/dx) - (test_old_x[b]/dx) == 0) and
                        (test_old_y[b+1]/dy) - (test_old_y[b]/dy) == -1):
                            
                            length_old += 1
                    
                    elif (((test_old_x[b+1]/dx) - (test_old_x[b]/dx) == 0)
                        and (test_old_y[b+1]/dy) - (test_old_y[b]/dy) == 1):
                    
                            length_old += 1
                            
                    elif (((test_old_x[b+1]/dx) - (test_old_x[b]/dx) == 1)
                        and (test_old_y[b+1]/dy) - (test_old_y[b]/dy) == 0):

                            length_old += 1
                    
                    elif (((test_old_x[b+1]/dx) - (test_old_x[b]/dx) == 1)
                        and (test_old_y[b+1]/dy) - (test_old_y[b]/dy) == -1):
                            
                            length_old += math.sqrt(2)

                    elif (((test_old_x[b+1]/dx) - (test_old_x[b]/dx) == 1)
                        and (test_old_y[b+1]/dy) - (test_old_y[b]/dy) == 1):
                            
                            length_old += math.sqrt(2)

                # if new river course < length of old
                # river course, then an avulsion will occur
                length_new_sum = sum(length_new)
                if sum(length_new) < length_old:
                    
                    loc = [a]         # avulsion location
                    avulsion_type = 1 # sets avulsion to be regional, may be 
                                        # updated again below (if local)
                
                    # maybe this should be len(test_old_x)-1?
                    for d in range(1,len(test_old_x)):
                        
                        x_diff = new_riv_x[-1] - riv_x[a+d]
                        y_diff = new_riv_y[-1] - riv_y[a+d]
                        
                        if x_diff == 0 and y_diff == 0:
                            
                            avulsion_type = 2   # local avulsion
                            
                            riv_x = new_riv_x + riv_x[a+d+1:]
                            riv_y = new_riv_y + riv_y[a+d+1:]
                            """
                            above doesn't change river mouth location unless it's
                            a regional avulsion
                            """ 

                            break
                    
                    if avulsion_type == 1: 
                    
                        riv_x = new_riv_x
                        riv_y = new_riv_y

                        n = downcut.cut_new(dx, dy, riv_x, riv_y, n, length_new,
                                        current_SL, a, ch_depth)

                    return (riv_x, riv_y, loc, SEL, SER, n, dn_fp, avulsion_type,
                            length_new_sum, length_old
                            )

                else:

                    if splay_type > 0:
                        n, dn_fp = FP.dep_splay(dy, dx, imax, jmax,
                                   riv_x, riv_y, new_riv_x, new_riv_y,
                                   ch_depth, n, a, dn_fp, splay_type,
                                   splay_dep)

            # if shortest path is not an avulsion criterion, then the new
            # steepest descent path will become the new course regardless
            # of new course length relative to the old course
            if short_path == 0:

                riv_x = new_riv_x
                riv_y = new_riv_y
                loc = [a]

    return (riv_x, riv_y, loc, SEL, SER, n, dn_fp, avulsion_type, length_new_sum,
            length_old)
