#!/usr/bin/env python
"""
author: Markus Stricker

Generic functions for averaging datasets.
"""

import numpy as np


def av_along_given_dir(data,ax):
    """
    Returns the average of an array along a given direction. A check is
    performed if the data is indexed in three dimensions. Assuming that
    the first three components are indices of voxels, the rest is the
    rest is the stuff which has to be averaged.

    Args:
        data (np.array): Input data of any size at lease larger than
            one entry per dimension.
        ax (int): Indicate axis for which to average.
    """

    if ax < 0 and ax > 3:
        print 'averaging.av_along_given_dir(): Specified wrong direction,\
               exiting.'
        exit

    return np.mean(data, axis=ax)

def radial_average(data):
    """
    Calculates the radial average of a two dimensional array.

    Args:
        data (np.array): Array with values, should be N==M, square cross
            section

    Returns:
        radial_avrg (np.array): Aray with number of points for averaging
            in first column and actual values in second column.
    """

    
    n = data.shape[0]
    cntr = np.array([n/2, n/2])
    radial_avrg = np.zeros([n/2,2])

    for i in range(n):
        for j in range(n):
            
            point = np.array([i,j])
            
            dist = np.linalg.norm(point-cntr)

            idx = int(dist)

            if idx >= n/2:
                continue
            
            radial_avrg[idx,1] += data[i][j]
            radial_avrg[idx,0] += 1.
    
    for i in range(n/2):
        radial_avrg[i,1] /= float(radial_avrg[i,0])

    return radial_avrg
            
    
def movingaverage(interval, window_size):
    """
    Returns a smoothed version of the input data.
    This is done by convoluting the data with a window.
    
    Args:
        interval (np.array): Array with values
        window_size (int) : Size of the averaging window, number
            of datapoints used

    Returns:
        conv (np.array): Same array, but smooth
    """

	
    window = np.ones(int(window_size))/float(window_size)
    bound_corr = interval[-window_size:]
    conv = np.convolve(interval, window, 'same')
    conv[-window_size:] = bound_corr

    return conv


def mirror_symmetry_average(data):
    """
    Returns a symmetric averga (as e.g. in crystal symmetry).

    Args:
        data (np.array NX1): Array with values

    Returns:
        mirror averaged value
    """

    return (data+data[::-1])/2.
