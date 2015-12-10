#!/usr/bin/env python
"""
author: Markus Stricker

Generic functions for averaging datasets.
"""

import numpy as np


def av_along_given_dir(data,direction):
    """
    Returns the average of an array along a given direction. A check is
    performed if the data is indexed in three dimensions.

    Args:
        data (np.array): Input data of any size at lease larger than
            one entry per dimension.
        direction (string): Indicate direction with string 'x', 'y' or 'z' 
    """

    if direction=='x':
        idx = 0
    elif direction=='y':
        idx = 1
    elif direction=='z':
        idx = 2
    else:
        print 'This direction does not exist:', direction
        return

    # Sanity checks for size

    

