#! /usr/local/bin/python
import numpy as np


def make_profile(n, riv_i, riv_j):
    """Get elevations along the river profile."""
    return n[riv_i, riv_j]
