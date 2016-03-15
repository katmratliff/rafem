#! /usr/local/bin/python
import yaml
import numpy as np
import pdb

def read_params_from_file(fname):
    """Read model parameters from a file.

    Parameters
    ----------
    fname : str
        Name of YAML-formatted parameters file.

    Returns
    -------
    dict
        A dict of parameters for the heat model.
    """
    with open(fname, 'r') as fp:
        params = yaml.load(fp)

#    params.setdefault('shape', (10, 20))
#    params.setdefault('spacing', (1., 1.))
#    params.setdefault('alpha', 1)
#    params.setdefault('end_time', 100.)

    return params


def is_diagonal_neighbor(sub0, sub1):
    return (sub0[0] != sub1[0]) and (sub0[1] != sub1[1])


def is_same_row(sub0, sub1):
    return sub0[0] == sub1[0]


def channel_is_superelevated(z, sub, channel_depth, super_ratio):
    """Check if a channel location is super-elevated.

    Parameters
    ----------
    z : ndarray
        2D-array of elevations.
    sub : tuple of int
        Row-column subscripts into *z*.
    channel_depth : float
        The depth of the channel.
    super_ratio : float
        Ratio of bankfull elevation to out-of-channel elevation for
        channel to be super-elevated.

    Returns
    -------
    boolean
        `True` if channel is super-elevated. Otherwise, `False`.
    """
    z_bankfull = z[sub] + channel_depth
    z_left, z_right = z[sub[0], sub[1] - 1], z[sub[0], sub[1] + 1]

    threshold = super_ratio * channel_depth

    return (z_bankfull - z_right >= threshold or
            z_bankfull - z_left >= threshold)


def get_link_lengths(path, dx=1., dy=1.):
    from itertools import izip

    DIAGONAL_LENGTH = np.sqrt(dx ** 2. + dy ** 2.)

    lengths = np.empty(len(path[0]) - 1, dtype=np.float)
    ij_last = path[0][0], path[1][0]
    for n, ij in enumerate(izip(path[0][1:], path[1][1:])):
        if is_diagonal_neighbor(ij, ij_last):
            lengths[n] = DIAGONAL_LENGTH
        elif is_same_row(ij, ij_last):
            lengths[n] = dx
        else:
            lengths[n] = dy
        ij_last = ij
    return lengths


def get_channel_distance(path, dx=1., dy=1.):
    total_distance = get_link_lengths(path, dx=dx, dy=dy).cumsum()
    return np.append(0, total_distance)


def find_path_length(n, path, sea_level, ch_depth, slope, dx=1., dy=1.):
    beach_len = find_beach_length(n, (path[0][-2], path[1][-2]),
                                  (path[0][-1], path[1][-1]), sea_level,
                                  ch_depth, slope, dx=dx, dy=dy)

    lengths = (get_link_lengths(path, dx=dx, dy=dy))
    lengths[-1] = np.divide(lengths[-1], 2.) + beach_len
    riv_length = lengths.sum()

    return (riv_length)


def find_point_in_path(path, sub):
    """Find a point in a path.

    Parameters
    ----------
    path : tuple of array_like
        Tuple of arrays of int that represent indices into a matrix.
    sub : tuple in int
        Indices to search *path* for.

    Returns
    -------
    int or None
        Index into *path* if the indices are found. Otherwise, ``None``.

    Examples
    --------
    >>> path = ((0, 1, 4, 6, 5), (3, 1, 7, 8, 10))
    >>> find_point_in_path(path, (4, 7))
    2
    >>> find_point_in_path(path, (99, 2)) is None
    True
    """
    try:
        return zip(*path).index(sub)
    except ValueError:
        return None


def set_linear_profile(z, riv_ij, dx=1., dy=1.):
    """Set elevations along a path to be linear.

    Examples
    --------
    >>> import numpy as np
    >>> z = np.array([[0, 0, 0, 0, 0],
    ...               [1, 2, 1, 1, 1],
    ...               [2, 2, 2, 2, 2],
    ...               [3, 3, 3, 3, 3]], dtype=float)
    >>> ij = [(0, 2), (1, 2), (2, 2), (3, 2)]
    >>> set_linear_profile(z, ij)
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 1.,  2.,  1.,  1.,  1.],
           [ 2.,  2.,  2.,  2.,  2.],
           [ 3.,  3.,  3.,  3.,  3.]])

    >>> ij = [(0, 0), (1, 1), (2, 2), (3, 1)]
    >>> set_linear_profile(z, ij, dx=3., dy=4.)
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 1.,  1.,  1.,  1.,  1.],
           [ 2.,  2.,  2.,  2.,  2.],
           [ 3.,  3.,  3.,  3.,  3.]])

    """
    z0 = z[riv_ij[0]]
    dz = z[riv_ij[-1]] - z[riv_ij[0]]

    lengths = get_link_lengths(zip(*riv_ij), dx=dx, dy=dy)
    ds = lengths.sum()

    # z[zip(*riv_ij[1:])] = z0 + dz / ds * lengths.cumsum()

    z[zip(*riv_ij[:-1])] = z[riv_ij[-1]] + np.arange(len(riv_ij) - 1, 0, -1) * 1e-6

    return z


def fill_upstream(z, riv_ij, dx=1., dy=1.):
    """Fill depressions upstream of a pit.

    Exammples
    ---------
    >>> import numpy as np
    >>> z = np.array([[3, 3, 3, 3, 3],
    ...               [2, 2, 2, 2, 2],
    ...               [1, 1, 0, 1, 1],
    ...               [0, 0, 1, 0, 0]], dtype=float)
    >>> ij = [(0, 2), (1, 2), (2, 2), (3, 2)]
    >>> fill_upstream(z, ij)
    array([[ 3. ,  3. ,  3. ,  3. ,  3. ],
           [ 2. ,  2. ,  2. ,  2. ,  2. ],
           [ 1. ,  1. ,  1.5,  1. ,  1. ],
           [ 0. ,  0. ,  1. ,  0. ,  0. ]])

    >>> z = np.array([[3, 3, 1, 3, 3],
    ...               [2, 2, 0, 2, 2],
    ...               [1, 1, 1, 1, 1],
    ...               [0, 0, 2, 0, 0]], dtype=float)
    >>> ij = [(0, 2), (1, 2), (2, 2), (3, 2)]
    >>> fill_upstream(z, ij)
    """
    fill_to = z[riv_ij[-1]]

    for n in xrange(len(riv_ij) - 2, -1, -1):
        if z[riv_ij[n]] > fill_to:
            break
        # z[riv_ij[n]] = z[riv_ij[n + 1]] + 1e-6

    if z[riv_ij[n]] <= fill_to:
        z[riv_ij[n]] = fill_to + 1e-6

    set_linear_profile(z, riv_ij[n:], dx=dx, dy=dy)

    return z


def sort_lowest_neighbors(n, sub):
    """Sort neighbor cells by elevation."""
    i, j = sub

    if j == n.shape[1] - 1:
        di, dj  = np.array([0, 1, 1]), np.array([-1, -1, 0])
    elif j == 0:
        di, dj  = np.array([1, 1, 0]), np.array([0, 1, 1])
    else:
        di, dj = np.array([0, 1, 1, 1, 0]),  np.array([-1, -1, 0, 1, 1])

    sorted = np.argsort(n[i + di, j + dj])

    return i + di[sorted], j + dj[sorted]

def find_beach_length(n, sub0, sub1, sea_level, channel_depth, slope, dx=1., dy=1.):
    """Find length of beach in shoreline cell.

    Parameters
    ----------
    z : ndarray
        2D-array of elevations.
    sub0 : tuple of int
        Row-column subscripts into *z*. (cell prior to beach)
    sub1 : tuple of int
        Row-column subscripts into *z*. (beach cell)
    sea_level : float
        Elevation of current sea level.
    channel_depth : float
        The depth of the channel.

    Returns
    -------
    Length of beach in a cell. """

    cell_elev = n[sub1] + channel_depth - sea_level

    max_cell_h = slope * dx

    DIAGONAL_LENGTH = np.sqrt(dx ** 2. + dy ** 2.)
    
    if is_diagonal_neighbor(sub0, sub1):
        d_dist = DIAGONAL_LENGTH

    elif is_same_row(sub0, sub1):
        d_dist = dx

    else:
        d_dist = dy

    return (cell_elev / max_cell_h) * d_dist


def lowest_cell_elev(n, sub):
    i,j = sub

    if j == 0 and i == 0:
        di, dj = np.array([1, 1, 0]), np.array([0, 1, 1])
    elif j == 0 and i == n.shape[0] - 1:
        di, dj = np.array([-1, -1, 0]), np.array([0, 1, 1])
    elif j == n.shape[1] - 1 and i == 0:
        di, dj = np.array([0, 1, 1]), np.array([-1, -1, 0])
    elif j == n.shape[1] - 1 and i == n.shape[0] - 1:
        di, dj = np.array([0, -1, -1]), np.array([-1, -1, 0])
    elif j == n.shape[1] - 1:
        di, dj  = np.array([-1, -1, 0, 1, 1]), np.array([0, -1, -1, -1, 0])
    elif j == 0:
        di, dj  = np.array([-1, -1, 0, 1, 1]), np.array([0, 1, 1, 1, 0])
    elif i == n.shape[0] - 1:
        di, dj = np.array([0, -1, -1, -1, 0]), np.array([-1, -1, 0, 1, 1])
    elif i == 0:
        di, dj = np.array([0, 1, 1, 1, 0]), np.array([-1, -1, 0, 1, 1])
    else:
        di, dj = np.array([0, -1, -1, -1, 0, 1, 1, 1]),  np.array([-1, -1, 0, 1, 1, 1, 0, -1])

    lowest = np.amin(n[i + di, j + dj])

    return lowest


def lowest_face(n, sub):
    i,j = sub

    shore_cell = 0

    if j == 0 and i == 0:
        di, dj = np.array([1, 0]), np.array([0, 1])
    elif j == 0 and i == n.shape[0] - 1:
        di, dj = np.array([-1, 0]), np.array([0, 1])
    elif j == n.shape[1] - 1 and i == 0:
        di, dj = np.array([0, 1]), np.array([-1, 0])
    elif j == n.shape[1] - 1 and i == n.shape[0] - 1:
        di, dj = np.array([0, -1]), np.array([-1, 0])
    elif j == n.shape[1] - 1:
        di, dj  = np.array([-1, 0, 1]), np.array([0, -1, 0])
    elif j == 0:
        di, dj  = np.array([-1, 0, 1]), np.array([0, 1, 0])
    elif i == n.shape[0] - 1:
        di, dj = np.array([0, -1, 0]), np.array([-1, 0, 1])
    elif i == 0:
        di, dj = np.array([0, 1, 0]), np.array([-1, 0, 1])
    else:
        di, dj = np.array([0, -1, 0, 1]),  np.array([-1, 0, 1, 0])

    lowest_face = np.amin(n[i + di, j + dj])

    return lowest_face


def fix_elevations(z, riv_i, riv_j, ch_depth, sea_level, slope, dx, max_rand):

    test_elev = z - sea_level
    max_cell_h = slope * dx
    riv_prof = test_elev[riv_i, riv_j]
    test_elev[riv_i, riv_j] += ch_depth

    riv_cells = np.zeros_like(z)
    riv_cells[riv_i, riv_j] = 1

    for i in xrange(1, test_elev.shape[0]):
        for j in xrange(test_elev.shape[1]):
            if riv_cells[i,j]:
                break
            # fix lakes (why do they form????) and partially full cells that aren't
            # part of the shoreline
            if (test_elev[i,j] < max_cell_h and
                np.any(test_elev[:i,j]) > 0 and
                np.any(test_elev[i+1:,j]) > 0 and
                np.any(test_elev[i,:j]) > 0 and
                np.any(test_elev[i,j+1:] > 0)):

                if test_elev[i,j] <= 0:
                    if j == 0:
                        test_elev[i,j] = test_elev[i,j+1] - np.random.rand()*max_rand
                    elif j == z.shape[1] - 1:
                        test_elev[i,j] = test_elev[i,j-1] - np.random.rand()*max_rand
                    else:
                        test_elev[i,j] = ((test_elev[i,j+1] + test_elev[i,j-1])/2
                                          - np.random.rand()*max_rand)
                else: test_elev[i,j] = max_cell_h + np.random.rand()*max_rand
            # Note: below keeps things sloping seaward. Needs revision to slope 
            # towards nearest shoreline cell (or perhaps a better idea is out there?)
            if riv_cells[i-1,j]:
                break
            if test_elev[i,j] > 0:
                if test_elev[i,j] >= test_elev[i-1,j] and lowest_face(test_elev,(i,j)) > 0:
                    test_elev[i-1,j] = test_elev[i,j] + (np.random.rand() * (slope*0.1))
    
    test_elev[riv_i, riv_j] = riv_prof

    z = test_elev + sea_level

    return z

    # for i in xrange(1,test_elev.shape[0]):
    #     for j in xrange(test_elev.shape[1]):
    #         if riv_cells[i,j]:
    #             break
    #         if test_elev[i,j] == 0:
    #             test_elev[i,j] == np.random.rand() * (slope*0.1)
    #         if test_elev[i,j] < 0 and lowest_cell_elev(test_elev, (i,j)) >= 0:
    #             #pdb.set_trace()
    #             test_elev[i,j] = slope*dx + np.random.rand()*max_rand
    #         # NEED TO HAVE SOMETHING HERE THAT WILL ROUTE WATER TO SEA??
    #         # perhaps find potential shoreline configurations and fix that way
    #         if riv_cells[i-1,j]:
    #                 break
    #         if test_elev[i,j] >= test_elev[i-1,j]:
    #             if test_elev[i,j] >= 0 and test_elev[i-1,j] >= 0:
    #                 #pdb.set_trace()
    #                 test_elev[i-1,j] = test_elev[i,j] + (np.random.rand() * (slope*0.1))

    # test_elev[riv_i, riv_j] = riv_prof

    # z = test_elev + sea_level

    # return z


