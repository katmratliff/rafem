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

    params.setdefault('shape', (10, 20))
    params.setdefault('spacing', (1., 1.))
    params.setdefault('alpha', 1)
    params.setdefault('end_time', 100.)

    return params