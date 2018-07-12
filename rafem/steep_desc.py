#! /usr/local/bin/python
import warnings

import numpy as np
from six.moves import range

from . import downcut
from .avulsion_utils import lowest_cell_elev, sort_lowest_neighbors


def lowest_neighbor(n, sub):
    """Find lowest neighbor value around a point.
    Parameters
    ----------
    n : ndarray
        Grid of elevations.
    sub : tuple of int
        Row/column subscript into elevation matrix.
    Returns
    -------
    tuple of int
        Subscripts into elevation matrix of the lowest neighbor.
    """
    i, j = sub

    if j == n.shape[1] - 1:
        di, dj  = np.array([0, 1, 1]), np.array([-1, -1, 0])
    elif j == 0:
        di, dj  = np.array([1, 1, 0]), np.array([0, 1, 1])
    else:
        di, dj = np.array([0, 1, 1, 1, 0]),  np.array([-1, -1, 0, 1, 1])

    lowest = np.argmin(n[i + di, j + dj])

    return i + di[lowest], j + dj[lowest]


def lowest_neighbor_prograde(n, sub):
    i, j = sub

    if j == n.shape[1] - 1:
        di, dj  = np.array([0, 1, 1]), np.array([-1, -1, 0])
    elif j == 0:
        di, dj  = np.array([1, 1, 0]), np.array([0, 1, 1])
    else:
        di, dj = np.array([0, 1, 1, 1, 0]),  np.array([-1, -1, 0, 1, 1])

    subaerial_cells = np.where(n[i + di, j + dj] > 0)

    subaerial_low = min(x for x in (n[i + di, j + dj]) if x > 0)
    lowest = np.where((n[i + di, j + dj]) == subaerial_low)[0][0]

    return i + di[lowest], j + dj[lowest]


def below_sea_level(z, sea_level):
    """Check if an elevation is below sea level.
    Parameters
    ----------
    z : float
        Elevation.
    sea_level : float
        Elevation of sea level.
    Returns
    -------
    boolean
        `True` if at or below sea level. Otherwise, `False`.
    """
    return z <= sea_level


def at_river_mouth(z, sub, z0):
    """Check if a cell is at the river mouth.
    Parameters
    ----------
    z : ndarray
        2D-array of elevations.
    sub : tuple of int
        Row and column subscript into *z*.
    z0 : float
        Elevation of sea level (or `None`).
    Returns
    -------
    boolean
        True if the cell at the given subscript is at the river mouth.
    """
    try:
        return sub[0] == z.shape[0] - 1 or below_sea_level(z[sub], z0)
    except IndexError:
        return True

def at_end_of_domain(z, sub):
    """Check if a cell a river mouth at the end of domain.
    Parameters
    ----------
    z : ndarray
        2D-array of elevations.
    sub : tuple of int
        Row and column subscript into *z*.
    Returns
    -------
    boolean
        True if the cell at the given subscript is at the river mouth.
    """
    try:
        return sub[0] == z.shape[0] - 1
    except IndexError:
        return True

def riv_cell_at_sea_level(z, sub, z0):
    """Check if a river cell is at sea level.
    Parameters
    ----------
    z : ndarray
        2D-array of elevations.
    sub : tuple of int
        Row and column subscript into *z*.
    z0 : float
        Elevation of sea level (or `None`).
    Returns
    -------
    boolean
        True if the cell at the given subscript is at the river mouth.
    """
    try:
        below_sea_level(z[sub], z0)
    except IndexError:
        return True

def find_course(z, riv_i, riv_j, SE_loc, channel_depth, sea_level=None):
    """Find the course of a river.
    Given a river course as subscripts, (*i*, *j*), into an array of
    elevations, (*z*), find the river path until the coast (*sea_level*) or
    the end of the domain.
    Parameters
    ----------
    z : ndarray
        Grid elevations.
    riv_i : array_like of int
        Row subscripts (into *n*) for river.
    riv_j : array_like of int
        Column subscripts (into *n*) for river.
    sea_level : None, optional
        The current sea level.
    Returns
    -------
    tuple of array_like
        Row and column indices of the new river path.
    Examples
    --------
    >>> import numpy as np
    >>> z = np.array([[4., 3., 4.],
    ...               [2., 3., 3.],
    ...               [2., 1., 2.]])
    >>> riv_i, riv_j = np.zeros(9, dtype=int), np.zeros(9, dtype=int)
    >>> riv_i[0], riv_j[0] = 0, 1
    >>> find_course(z, riv_i[:1], riv_j[:1], 1)
    (array([0, 1, 2]), array([1, 0, 1]))
    >>> find_course(z, riv_i[:1], riv_j[:1], sea_level=2.)
    (array([0, 1]), array([1, 0]))
    >>> z = np.array([[4., 3., 4.],
    ...               [2., 3., 3.],
    ...               [2., 1., 2.],
    ...               [2., 1.5, 2]])
    >>> find_course(z, riv_i[:1], riv_j[:1], sea_level=0.)
    (array([0, 1, 2, 3]), array([1, 0, 1, 1]))
    >>> z
    """
    # function to find the steepest descent route
    # note: this needs to be improved to remove potential bias that may occur
    # if two or more cells have the steepest descent elevation
    old_course = zip(riv_i, riv_j)

    n_levee = np.copy(z)
    n_levee[riv_i, riv_j] += channel_depth

    riv_i = riv_i[:SE_loc]
    riv_j = riv_j[:SE_loc]

    assert(riv_i.size > 0 and riv_j.size > 0)

    if sea_level is None:
        sea_level = - np.finfo(float).max

    for n in range(1, riv_i.size):
        if at_end_of_domain(z, (riv_i[n - 1], riv_j[n - 1])):
            return riv_i[:n], riv_j[:n]

    # for n in range(1, riv_i.size):
    #     if riv_cell_at_sea_level(z, (riv_i[n - 1], riv_j[n - 1]), sea_level):
    #         return riv_i[:n-1], riv_j[:n-1]

    new_i = np.empty(z.size, dtype=np.int)
    new_j = np.empty(z.size, dtype=np.int)

    new_i[:len(riv_j)] = riv_i[:]
    new_j[:len(riv_i)] = riv_j[:]

    pits = True
    while pits:
        for n in range(riv_i.size, new_i.size):
            # if at_river_mouth(z, (new_i[n - 1], new_j[n - 1]), sea_level):
            #     pits = False
            #     break

            if at_end_of_domain(z, (new_i[n - 1], new_j[n - 1])):
                pits = False
                break

            sorted_n = sort_lowest_neighbors(n_levee, (new_i[n - 1], new_j[n - 1]))

            if (sorted_n[0][0], sorted_n[1][0]) not in zip(new_i[:n - 1], new_j[:n - 1]):
                downstream_ij = (sorted_n[0][0], sorted_n[1][0])
            elif (sorted_n[0][1], sorted_n[1][1]) not in zip(new_i[:n - 1], new_j[:n - 1]):
                downstream_ij = (sorted_n[0][1], sorted_n[1][1])
            elif (sorted_n[0][2], sorted_n[1][2]) not in zip(new_i[:n - 1], new_j[:n - 1]):
                downstream_ij = (sorted_n[0][2], sorted_n[1][2])
            else:
                raise RuntimeError('river course is going crazy!')

            if downstream_ij not in old_course and below_sea_level(z[downstream_ij], sea_level):
                pits = False
                break

            if downstream_ij in old_course:
                new_i[n], new_j[n] = downstream_ij
                n += 1
                pits = False
                break

            if z[downstream_ij] > z[new_i[n - 1], new_j[n - 1]]:
                new_i[n], new_j[n] = downstream_ij
                z[new_i[n - 1], new_j[n - 1]] +=  1e-6
            else:
                new_i[n], new_j[n] = downstream_ij

            # new_i[n], new_j[n] = lowest_neighbor(z, (new_i[n - 1], new_j[n - 1]))

    if n == 0:
        raise RuntimeError('new river length is zero!')

    return new_i[:n], new_j[:n]


def update_course(z, riv_i, riv_j, ch_depth, slope, sea_level=None, dx=1., dy=1.):

    if sea_level is None:
        sea_level = - np.finfo(float).max

    course_update = 0

    last_elev = z[riv_i[-1], riv_j[-1]] + ch_depth - sea_level
    max_cell_h = slope * dx

    test_elev = np.copy(z)
    test_elev -= sea_level
    test_elev[riv_i, riv_j] += 2 * ch_depth

    low_adj_cell = lowest_cell_elev(test_elev, (riv_i[-1], riv_j[-1]))

    # check for coastal avulsion (happens if river is prograding too far alongshore)
    if ((riv_i[-1] == riv_i[-2] == riv_i[-3] == riv_i[-4] == riv_i[-5]) and   # if last 5 river cells flowing alongshore
        # (z[riv_i[-1+1],riv_j[-1]] > sea_level) and  # if land between river and ocean for last 3 cells
        (z[riv_i[-2]+1,riv_j[-2]] > sea_level) and
        (z[riv_i[-3]+1,riv_j[-3]] > sea_level) and
        (z[riv_i[-4]+1,riv_j[-4]] > sea_level) and
        (z[riv_i[-5]+1,riv_j[-5]] > sea_level) and
        (z[riv_i[-1]+2,riv_j[-1]] <= sea_level) and
        (z[riv_i[-2]+2,riv_j[-2]] <= sea_level) and
        (z[riv_i[-3]+2,riv_j[-3]] <= sea_level) and
        (z[riv_i[-4]+2,riv_j[-4]] <= sea_level) and
        (z[riv_i[-5]+2,riv_j[-5]] <= sea_level)):

        # fill up old river mouth
        z[riv_i[-1],riv_j[-1]] += ch_depth

        # turn river towards ocean
        riv_i[-1] = riv_i[-2]+1
        riv_j[-1] = riv_j[-2]

        if (z[riv_i[-1], riv_j[-1]] - sea_level) < (0.001 * max_cell_h):
            z[riv_i[-1], riv_j[-1]] = (0.001 * max_cell_h) + sea_level

        # z[riv_i[-1],riv_j[-1]] -= ch_depth

        downcut.cut_new(riv_i[-3:], riv_j[-3:], z, sea_level, ch_depth,
                        dx=dx, dy=dy)

        course_update = 7   # coastal avulsion

    elif last_elev <= 0:
        riv_i = riv_i[:-1]
        riv_j = riv_j[:-1]
        course_update = 4   # shortened course

    # if river mouth surrounded by land
    elif low_adj_cell > 0:
        new_riv_i, new_riv_j = find_course(z, riv_i, riv_j, len(riv_i), 
                                           ch_depth, sea_level=sea_level)

        new_riv_length = new_riv_i.size - riv_i.size

        if new_riv_length > 0:
            riv_i = new_riv_i
            riv_j = new_riv_j

            if (z[riv_i[-1], riv_j[-1]] - sea_level) < (0.001 * max_cell_h):
                z[riv_i[-1], riv_j[-1]] = (0.001 * max_cell_h) + sea_level

            downcut.cut_new(riv_i[-(new_riv_length+2):], riv_j[-(new_riv_length+2):],
                                z, sea_level, ch_depth, dx=dx, dy=dy)

            course_update = 6 # lengthened land-locked course
        
        else:
            riv_i = riv_i
            riv_j = riv_j

    # if river mouth needs to prograde
    elif last_elev >= max_cell_h:       
        sorted_n = sort_lowest_neighbors(test_elev, (riv_i[-1], riv_j[-1]))
        subaerial_loc = np.where(test_elev[sorted_n] > 0)

        if len(subaerial_loc[0]):
            subaerial_cells = sorted_n[0][subaerial_loc], sorted_n[1][subaerial_loc]

            if (subaerial_cells[0][0], subaerial_cells[1][0]) not in zip(riv_i, riv_j):
                riv_i = np.append(riv_i, subaerial_cells[0][0])
                riv_j = np.append(riv_j, subaerial_cells[1][0])

                if (z[riv_i[-1], riv_j[-1]] - sea_level) < (0.001 * max_cell_h):
                    z[riv_i[-1], riv_j[-1]] = (0.001 * max_cell_h) + sea_level
                
                ### line below not needed if downcutting ###
                # z[riv_i[-1], riv_j[-1]] -= ch_depth

                downcut.cut_new(riv_i[-3:], riv_j[-3:], z, sea_level, ch_depth,
                                dx=dx, dy=dy)

                course_update = 5

            else:
                riv_i = riv_i
                riv_j = riv_j

        else:
            riv_i = riv_i
            riv_j = riv_j

    else:
        riv_i = riv_i
        riv_j = riv_j

    return riv_i, riv_j, course_update
