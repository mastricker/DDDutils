#!/usr/bin/env python

"""
Helper functions for sortin DDD data
"""


import numpy as np



def correctformonotony(x):
    delta = 0.

    for i in np.arange(x.shape[0]-1):
        delta_tmp = x[i+1] - x[i]
	x[i] += delta
        
        if (delta_tmp < 0. and np.absolute(delta_tmp)>0.0003):
	    delta += np.absolute(delta_tmp)			
            
    x[-1] += delta
            
    return x
