#! /usr/local/bin/python

import steep_desc
import downcut
import FP
import numpy as np
import math

from avulsion_utils import (find_point_in_path, channel_is_superelevated,
                            find_path_length)


def avulse_to_new_path(z, old, new, sea_level, channel_depth, dx=1., dy=1.):
    old_i, old_j = old
    new_i, new_j = new
    # sets avulsion to be regional, may be updated again below (if local)
            
    # maybe this should be len(test_old_x)-1?
    ind = find_point_in_path((old_i[1:], old_j[1:]), (new_i[-1], new_j[-1]))

    if ind is not None:
        avulsion_type = 2

        new_i = np.append(new_i, old_i[ind + 1:])
        new_j = np.append(new_j, old_j[ind + 1:])
    else:
        downcut.cut_new(new_i, new_j, z, sea_level, channel_depth, dx=dx,
                        dy=dy)

    return (new_i, new_j)


# determines if there is an avulsion along river course
def find_avulsion(riv_i, riv_j, n, super_ratio, current_SL, ch_depth,
                  short_path, splay_type, splay_dep, dx=1., dy=1.):
    new = riv_i, riv_j
    old = riv_i, riv_j

    for a in xrange(1, len(riv_i)):
        if channel_is_superelevated(n, (riv_i[a], riv_j[a]), ch_depth,
                                    super_ratio):

            # if superelevation greater than trigger ratio, determine
            # length of new steepest descent path

            new = steep_desc.find_course(n, riv_i[:a - 1], riv_j[:a - 1],
                                         sea_level=current_SL)

            # if using the shortest path as an avulsion criterion, then
            # the lengths of the previous and newly calculated paths will
            # be compared
            if short_path == 1:
                new_length = find_path_length(new, dx=dx, dy=dy)
                old_length = find_path_length(old, dx=dx, dy=dy)

                if new_length < old_length:
                    # if new river course < length of old
                    # river course, then an avulsion will occur
                    new = avulse_to_new_path(n, (riv_i[a - 1:], riv_j[a - 1:]),
                                             (new[0][a - 1:], new[1][a - 1:]),
                                             current_SL, ch_depth, dx=dx,
                                             dy=dy)

                    new = (np.append(riv_i[:a - 1], new[0]),
                           np.append(riv_j[:a - 1], new[1]))

                elif splay_type > 0:
                    FP.dep_splay(n, (new[0][a], new[1][a]), (riv_i, riv_j),
                                 splay_dep, splay_type=splay_type)
            # if shortest path is not an avulsion criterion, then the new
            # steepest descent path will become the new course regardless
            # of new course length relative to the old course

    return new
