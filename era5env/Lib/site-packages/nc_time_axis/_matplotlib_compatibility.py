from numbers import Number

import numpy as np
from numpy import ma


def is_numlike(x):
    """
    The Matplotlib datalim, autoscaling, locators etc work with scalars
    which are the units converted to floats given the current unit.  The
    converter may be passed these floats, or arrays of them, even when
    units are set.

    Vendored for matplotlib.  This function will not be needed once the
    minimum version of matplotlib supported by nc-time-axis is at least
    3.5.  See GitHub issue 97 for more details.
    """
    if np.iterable(x):
        for thisx in x:
            if thisx is ma.masked:
                continue
            return isinstance(thisx, Number)
    else:
        return isinstance(x, Number)
