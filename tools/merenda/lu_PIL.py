# Polarity Inversion Line computation
#
#
#
#
#
#
#

import sunpy as sp
import sunpy.map as sunpymap
import astropy.units as u
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt


def sign(bo):

    # Sign function: maps positive values to 1, 0 to 0 and negative values to -1
    if bo > 0:
        return 1
    elif bo == 0:
        return 0
    else:
        return -1


def line(x, y, a, b, c):

    return a * x + b * y + c
