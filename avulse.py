#! /usr/local/bin/python

import steep_desc
import downcut
import FP
import numpy as np
import math


def is_diagonal_neighbor(sub0, sub1):
    return sub0[0] != sub1[0] and sub0[1] != sub1[1]


# determines if there is an avulsion along river course
def find_avulsion(riv_i, riv_j, n, super_ratio, current_SL, ch_depth,
                  short_path, dn_fp, splay_type, splay_dep):

    loc = []
    SEL = np.zeros(len(riv_i))
    SER = np.zeros(len(riv_j))
    avulsion_type = 0
    length_new_sum = 0    
    length_old = 0

    for a in xrange(1, len(riv_i)):

        ch_Z = n[riv_i[a], riv_j[a]] + ch_depth   # bankfull elev.
        LHS = n[riv_i[a], riv_j[a] - 1]
        RHS = n[riv_i[a], riv_j[a] + 1]

        # normalized superelevation ratio on left side
        SEL[a] = ((ch_Z - LHS) / ch_depth)

        # normalized superelevation ratio on right side
        SER[a] = ((ch_Z - RHS) / ch_depth)

        if SEL[a] >= super_ratio or SER[a] >= super_ratio:

            # if superelevation greater than trigger ratio, determine
            # length of new steepest descent path
            new_riv_i = riv_i[:a-1].copy()
            new_riv_j = riv_j[:a-1].copy()

            steep_desc.find_course(n, new_riv_i, new_riv_j,
                                   sea_level=current_SL)

            # if using the shortest path as an avulsion criterion, then
            # the lengths of the previous and newly calculated paths will
            # be compared
            if short_path == 1:
                
                # duplicates arrays so that length can be compared below
                test_new_i = new_riv_i[a:]
                test_new_j = new_riv_j[a:]
                test_old_i = riv_i[a:]
                test_old_j = riv_j[a:]
                length_new = []
                
                for c in xrange(len(test_new_i) - 1):
                    ij_cur = test_new_i[c], test_new_j[c]
                    ij_next = test_new_i[c + 1], test_new_j[c + 1]

                    if is_diagonal_neighbor(ij_cur, ij_next):
                        length_new.append(np.sqrt(2.))
                    else:
                        length_new.append(1.)
                
                for b in xrange(len(test_old_i) - 1):
                    ij_cur = test_old_i[b], test_old_j[b]
                    ij_next = test_old_i[b + 1], test_old_j[b + 1]
                    
                    if is_diagonal_neighbor(ij_cur, ij_next):
                        length_old += np.sqrt(2.)
                    else:
                        length_old += 1.

                # if new river course < length of old
                # river course, then an avulsion will occur
                length_new_sum = sum(length_new)
                if sum(length_new) < length_old:
                    
                    loc = [a]         # avulsion location
                    avulsion_type = 1 # sets avulsion to be regional, may be 
                                      # updated again below (if local)
                
                    # maybe this should be len(test_old_x)-1?
                    for d in xrange(1, len(test_old_i)):
                        i_diff = new_riv_i[-1] - riv_i[a + d]
                        j_diff = new_riv_j[-1] - riv_j[a + d]
                        
                        if i_diff == 0 and j_diff == 0:
                            
                            avulsion_type = 2   # local avulsion
                            
                            riv_i = new_riv_i + riv_i[a + d + 1:]
                            riv_j = new_riv_j + riv_j[a + d + 1:]
                            # above doesn't change river mouth location unless it's
                            # a regional avulsion

                            break
                        
                        else: d += 1
                    
                    if avulsion_type == 1: 
                    
                        riv_i = new_riv_i
                        riv_j = new_riv_j

                        downcut.cut_new(riv_i, riv_j, n, np.array(length_new),
                                        current_SL, a, ch_depth)

                    return (riv_i, riv_j, loc, n, dn_fp, avulsion_type,
                            length_new_sum, length_old)

                else:
                    if splay_type > 0:
                        n, dn_fp = FP.dep_splay(riv_i, riv_j, new_riv_i,
                                                new_riv_j, ch_depth, n, a,
                                                dn_fp, splay_type, splay_dep)

            # if shortest path is not an avulsion criterion, then the new
            # steepest descent path will become the new course regardless
            # of new course length relative to the old course
            if short_path == 0:

                riv_i = new_riv_i
                riv_j = new_riv_j
                loc = [a]
        a += 1

    return (riv_i, riv_j, loc, n, dn_fp, avulsion_type, length_new_sum,
            length_old)
