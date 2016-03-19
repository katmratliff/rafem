#! /usr/local/bin/python

import steep_desc
import downcut
import FP
import numpy as np
import math
import pdb

from avulsion_utils import (find_point_in_path, channel_is_superelevated,
                            find_path_length)


def avulse_to_new_path(z, old, new, sea_level, channel_depth, avulsion_type,
                       slope, dx=1., dy=1.,):
    """Avulse the river to a new path.

    Given two river paths, *old* and *new*, avulse the river to a new river
    path. If the end point of the new path is contained in the old river
    path, the resulting path is the new path up until this point and then
    the old path. Otherwise, the resulting path is the new river path and
    will be downcut.

    Parameters
    ----------
    z : ndarray
        2D array of elevations.
    old : tuple of array_like
        Tuple of i and j indices (into *z*) for the old path.
    new : tuple of array_like
        Tuple of i and j indices (into *z*) for the new path.
    sea_level : float
        Elevation of sea level.
    channel_depth : float
        Depth of the channel.
    avulsion_type : {0, 1, 2, 3}
        The type of the avulsion.
    dx : float, optional
        Spacing of columns of *z*.
    dy : float, optional
        Spacing of rows of *z*.

    Returns
    -------
    tuple
        Tuple of the new river path (as i, j indices) and the, possibly
        changed, avulsion type.

    Examples
    --------
    The following example uses a grid that looks like::

        o  +  *  *
        *  o  +  *
        *  *  +  *
        *  *  o  *
        *  o  *  *
    
    The old path is marked by `o`, the new path but `+`. The paths overlap
    (2, 2).

    >>> import numpy as np
    >>> z = np.ones((5, 4), dtype=float)

    >>> old = np.array((0, 1, 2, 3, 4)), np.array((0, 1, 2, 2, 1))
    >>> new = np.array((0, 1, 2)), np.array((1, 2, 2))
    >>> (new, atype) = avulse_to_new_path(z, old, new, 0., 0., 0)

    The new path follows the new path until the common point and then
    follows the old path. The new avulsion type is now 2.

    >>> new
    (array([0, 1, 2, 3, 4]), array([1, 2, 2, 2, 1]))
    >>> atype
    2

    In this example the old and new paths do not overlap::

        o  +  *  *
        *  o  +  *
        *  *  o  +
        *  *  o  +
        *  o  *  +

    >>> old = np.array((0, 1, 2, 3, 4)), np.array((0, 1, 2, 2, 1))
    >>> new = np.array((0, 1, 2, 3, 4)), np.array((1, 2, 3, 3, 3))
    >>> (new, atype) = avulse_to_new_path(z, old, new, 0., 0., 0)

    The new path is now, in fact, the actual new path and the avulsion
    type is unchanged.

    >>> new
    (array([0, 1, 2, 3, 4]), array([1, 2, 3, 3, 3]))
    >>> atype
    0
    """
    old_i, old_j = old
    new_i, new_j = new
    # sets avulsion to be regional, may be updated again below (if local)
            
    # maybe this should be len(test_old_x)-1?
    ind = find_point_in_path((old_i, old_j), (new_i[-1], new_j[-1]))

    if ind is not None:
        avulsion_type = 2

        new_i = np.append(new_i, old_i[ind + 1:])
        new_j = np.append(new_j, old_j[ind + 1:])
    else:
        if (z[riv_i[-1], riv_j[-1]] - sea_level) < (0.001 * max_cell_h):
            z[riv_i[-1], riv_j[-1]] = (0.001 * max_cell_h) + sea_level
        
        downcut.cut_new(new_i, new_j, z, sea_level, channel_depth, slope,
                        dx=dx, dy=dy)

    return (new_i, new_j), avulsion_type


# determines if there is an avulsion along river course
def find_avulsion(riv_i, riv_j, n, super_ratio, current_SL, ch_depth,
                  short_path, splay_type, splay_dep, slope, dx=1., dy=1.):
    new = riv_i, riv_j
    old = riv_i, riv_j
    avulsion_type = 0
    a = 0

    for a in xrange(1, len(riv_i)-1):
        if channel_is_superelevated(n, (riv_i[a], riv_j[a]),
                                    (riv_i[a-1], riv_j[a-1]),
                                    ch_depth, super_ratio, current_SL):
            # if superelevation greater than trigger ratio, determine
            # length of new steepest descent path

            new = steep_desc.find_course(n, riv_i, riv_j, a, ch_depth,
                                         sea_level=current_SL)

            # if using the shortest path as an avulsion criterion, then
            # the lengths of the previous and newly calculated subaerial
            # paths will be compared
            if short_path == 1:
                new_length = find_path_length(n, new, current_SL, ch_depth,
                                              slope, dx=dx, dy=dy)
                old_length = find_path_length(n, old, current_SL, ch_depth,
                                              slope, dx=dx, dy=dy)

                if new_length < old_length:
                    # if new subaerial river course < length of old
                    # river course, then an avulsion will occur
                    avulsion_type = 1

                    new, avulsion_type = avulse_to_new_path(n,
                                             (riv_i[a - 1:], riv_j[a - 1:]),
                                             (new[0][a - 1:], new[1][a - 1:]),
                                             current_SL, ch_depth, avulsion_type,
                                             slope, dx=dx, dy=dy)

                    new = (np.append(riv_i[:a - 1], new[0]),
                           np.append(riv_j[:a - 1], new[1]))

                    # fill up old channel... could be some fraction in the future
                    # (determines whether channels are repellors or attractors)

                    new_profile = n[new[0], new[1]]
                    n[riv_i[a:], riv_j[a:]] += ch_depth
                    n[new[0], new[1]] = new_profile

                    break

                elif splay_type > 0:
                    avulsion_type = 3

                    # below should just be a??? not a-1???
                    FP.dep_splay(n, (new[0][a-1], new[1][a-1]), (riv_i, riv_j),
                                 splay_dep, splay_type=splay_type)
            # if shortest path is not an avulsion criterion, then the new
            # steepest descent path will become the new course regardless
            # of new course length relative to the old course

    return new, avulsion_type, a
