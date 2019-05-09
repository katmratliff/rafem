#! /usr/local/bin/python
import yaml
import numpy as np
from scipy.ndimage import measurements
from six.moves import range


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
        params = yaml.safe_load(fp)

    return params


def is_diagonal_neighbor(sub0, sub1):
    return (sub0[0] != sub1[0]) and (sub0[1] != sub1[1])


def is_same_row(sub0, sub1):
    return sub0[0] == sub1[0]


def channel_is_superelevated(z, riv, behind, channel_depth,
                             super_ratio, sea_level):
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
    superelev = 0
    z_bankfull = z[riv] + channel_depth

    # cross-shore river orientation
    if behind[0]+1 == riv[0]:
        if riv[1] == 0:
            adj1 = z[riv[0], riv[1] + 1]
            adj2 = z[riv[0], riv[1] + 1] # not sure if this needs to be included
        elif riv[1] == z.shape[1] - 1:
            adj1 = z[riv[0], riv[1] - 1]
            adj2 = z[riv[0], riv[1] - 1]
        else:
            adj1, adj2 = z[riv[0], riv[1] - 1], z[riv[0], riv[1] + 1]

    # alongshore river orientation
    if behind[0] == riv[0]:
        if riv[0] == z.shape[0] - 1:
            adj1 = z[riv[0] + 1, riv[1]]
            adj2 = z[riv[0] + 1, riv[1]]
        else:
            adj1, adj2 = z[riv[0] + 1, riv[1]], z[riv[0] - 1, riv[1]]

    threshold = super_ratio * channel_depth

    if (adj1 < sea_level) and (adj2 < sea_level):
        superelev = 0
    elif adj1 < sea_level:
        if z_bankfull - adj2 >= threshold:
            superelev = 1
    elif adj2 < sea_level:
        if z_bankfull - adj1 >= threshold:
            superelev = 1
    else:
        if ((z_bankfull - adj2) or (z_bankfull - adj1)) >= threshold:
            superelev = 1

    return superelev


def get_link_lengths(path, dx=1., dy=1.):
    from six.moves import zip

    DIAGONAL_LENGTH = np.sqrt(dx ** 2. + dy ** 2.)

    lengths = np.empty(len(path[0]) - 1, dtype=np.float)
    ij_last = path[0][0], path[1][0]
    for n, ij in enumerate(zip(path[0][1:], path[1][1:])):
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
    beach_len = find_new_beach_length(n, (path[0][-2], path[1][-2]),
                                  (path[0][-1], path[1][-1]), sea_level,
                                  dx=dx, dy=dy)

    lengths = (get_link_lengths(path, dx=dx, dy=dy))
    lengths[-1] = np.divide(lengths[-1], 2.) + beach_len
    riv_length = lengths.sum()

    return (riv_length)


def find_riv_path_length(n, path, sea_level, ch_depth, slope, dx=1., dy=1.):
    beach_len = find_beach_length_riv_cell (n, (path[0][-2], path[1][-2]),
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

    z[zip(*riv_ij[:-1])] = z[riv_ij[-1]] + np.arange(len(riv_ij) - 1, 0, -1) * 1e-6

    return z

def set_linear_slope(z, riv_ij, dx=1., dy=1.):
    """Set slope along a path to be linear."""

    if len(riv_ij) > 1:
        z0 = z[riv_ij[0]]
        dz = z[riv_ij[-1]] - z[riv_ij[0]]

        lengths = get_link_lengths(zip(*riv_ij), dx=dx, dy=dy)
        ds = lengths.sum()

        z[zip(*riv_ij[1:])] = z0 + dz / ds * lengths.cumsum()

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

    for n in range(len(riv_ij) - 2, -1, -1):
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

def find_new_beach_length(n, sub0, sub1, sea_level, dx=1., dy=1.):

    cell_elev = n[sub1] - sea_level

    DIAGONAL_LENGTH = np.sqrt(dx ** 2. + dy ** 2.)
    
    if is_diagonal_neighbor(sub0, sub1):
        d_dist = DIAGONAL_LENGTH

    elif is_same_row(sub0, sub1):
        d_dist = dx

    else:
        d_dist = dy

    return cell_elev * d_dist

def find_beach_length_riv_cell (n, sub0, sub1, sea_level, channel_depth,
                                slope, dx=1., dy=1.):
    """Find length of beach in shoreline cell with river.

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


def fix_elevations(z, riv_i, riv_j, ch_depth, sea_level, slope, dx, max_rand, SLRR):

    test_elev = z - sea_level
    max_cell_h = slope * dx
    riv_prof = test_elev[riv_i, riv_j]
    test_elev[riv_i, riv_j] += 2*ch_depth

    # set new subaerial cells to marsh elevation
    test_elev[test_elev == 0] = max_cell_h

    # make mask for depressions
    ocean_mask = test_elev < max_cell_h
    labeled_ponds, ocean = measurements.label(ocean_mask)

    # # fill in underwater spots that are below SL (?)
    # below_SL = [z <= sea_level]
    # underwater_cells, big_ocean = measurements.label(below_SL)
    # underwater_cells[underwater_cells == big_ocean] = 0
    # test_elev[underwater_cells > 0] = max_cell_h + SLRR + (np.random.rand() * max_rand)

    # create an ocean and shoreline mask
    ocean_and_shore = np.copy(labeled_ponds)

    # create an ocean mask
    # ocean_cells = np.copy(ocean_and_shore)
    # ocean_and_shore[test_elev > 0] = 0

    # create mask for pond cells and fix them
    area = measurements.sum(ocean_mask, labeled_ponds, index=np.arange(labeled_ponds.max() + 1))
    areaPonds = area[labeled_ponds]
    labeled_ponds[areaPonds == areaPonds.max()] = 0

    #finish creating ocean and shoreline mask
    ocean_and_shore[areaPonds != areaPonds.max()] = 0

    # something here to get rid of ocean cells
    test_elev[labeled_ponds > 0] = max_cell_h + SLRR + (np.random.rand() * max_rand)

    # raise cells close to sea level above it
    test_elev[(test_elev >= max_cell_h) & (test_elev <= (max_cell_h + SLRR))] = \
        (max_cell_h + SLRR + (np.random.rand() * max_rand))

    riv_buffer = np.zeros_like(test_elev)
    riv_buffer[riv_i, riv_j] = 1
    riv_buffer[riv_i[1:]-1, riv_j[1:]] = 1

    for i in range(1, test_elev.shape[0]-1):
        for j in range(test_elev.shape[1]):
            if (not ocean_and_shore[i, j]
                and not ocean_and_shore[i-1, j]
                and not ocean_and_shore[i+1, j]
                and not riv_buffer[i, j]):
                if test_elev[i+1, j] >= test_elev[i, j]:
                    test_elev[i, j] = test_elev[i+1, j] + (np.random.rand() * slope)
    
    test_elev[riv_i, riv_j] = riv_prof

    z = test_elev + sea_level

    return z


def fill_abandoned_channel(breach_loc, n, new, riv_i, riv_j, current_SL,
                           ch_depth, slope, dx):

    new_profile = n[new[0], new[1]]
    
    max_shorecell_h = current_SL + (slope * dx)
    riv_cells = np.zeros_like(n)
    riv_cells[riv_i,riv_j] = 1
    riv_cells[new[0],new[1]] = 1

    n[riv_i[-1],riv_j[-1]] += ch_depth

    if len(riv_i) - breach_loc > 1:
        for k in range(breach_loc, len(riv_i)-1):
            if (riv_j[k] == n.shape[1] - 1) or (riv_j[k] == 0):
                n[riv_i[k], riv_j[k]] = max_shorecell_h + (np.random.rand() * slope)
            elif (riv_cells[riv_i[k], riv_j[k]+1] and riv_cells[riv_i[k], riv_j[k]-1]):
                n[riv_i[k], riv_j[k]] = max_shorecell_h + (np.random.rand() * slope)
            elif riv_cells[riv_i[k], riv_j[k]+1]:
                if n[riv_i[k], riv_j[k]-1] <= max_shorecell_h:
                    n[riv_i[k], riv_j[k]] = max_shorecell_h + (np.random.rand() * slope)
                else:
                    n[riv_i[k], riv_j[k]] = (n[riv_i[k], riv_j[k]-1]
                                             + (np.random.rand() * slope))
            elif riv_cells[riv_i[k], riv_j[k]-1]:
                if n[riv_i[k], riv_j[k]+1] <= max_shorecell_h:
                    n[riv_i[k], riv_j[k]] = max_shorecell_h + (np.random.rand() * slope)
                else:
                    n[riv_i[k], riv_j[k]] = (n[riv_i[k], riv_j[k]+1]
                                             + (np.random.rand() * slope))
            elif ((n[riv_i[k], riv_j[k]+1] <= max_shorecell_h) or
                  (n[riv_i[k], riv_j[k]-1] <= max_shorecell_h)):
                  n[riv_i[k], riv_j[k]] = max_shorecell_h + (np.random.rand() * slope)
            else:
                n[riv_i[k], riv_j[k]] = (((n[riv_i[k], riv_j[k]+1]
                                           + n[riv_i[k], riv_j[k]-1])/2)
                                           + (np.random.rand() * slope))

    n[new[0], new[1]] = new_profile

