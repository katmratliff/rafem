import yaml
import numpy as np

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


def find_path_length(path, dx=1., dy=1.):
    return get_link_lengths(path, dx=dx, dy=dy).sum()


def find_point_in_path(path, sub):
    try:
        return zip(path).index(sub)
    except ValueError:
        return None
