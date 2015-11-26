#!/usr/bin/env python

"""
Helper functions for sorting DDD data
"""


import numpy as np


def checkmonotony(x):
    return all(x<y for x,y in zip(x, x[1:]))


def correctformonotony(x, sortidx=0):
    """
    Function for correcting data, which has been restarted.
    x : np.array
    index : otional index for sorting (int)
    
    Assuming the data is in the first index and the second index of x is the column to specify data.
    Example: if x = np.loadtxt('result_sigma_eps.dat'), you have this structure.
    """

    n = -1
    try:
        m,n = x.shape
    except Exception:
        # one dimensional array
        m = x.shape[0]

#    print 'm,n',m,n

    if n < 0:
        mono = checkmonotony(x)
        if mono:
            result = x
            return result
    else:
        mono = checkmonotony(x[:,sortidx])
        if mono:
            result = x
            return result

    
    print 'mono before?',mono

    # get first line of data
    if n > 0:
#        print 'case 1'
        x_corr = np.empty([1,n])
        x_corr[0,:] = x[0,:]
        ctr = 0
        #print x.shape, x_corr.shape,sortidx,m,x_corr
        for i in range(1,m):
            #print x[i,sortidx],x_corr[ctr,sortidx]
            if (x[i,sortidx] > x_corr[ctr,sortidx]) and np.mod(x[i,sortidx],10) == 0:
                ctr += 1
                x_corr = np.vstack([x_corr, x[i,:]])

        mono = checkmonotony(x_corr[:,0])
    else:
        x_corr = np.empty(1)
        x_corr[0] = x[0]
        ctr = 0
        for i in range(1,x.shape[0]):
            if x[i] > x_corr[ctr] and np.mod(x[i],10) == 0:
                ctr += 1
                x_corr = np.vstack([x_corr,x[i]])

    print 'array size before/after',x.shape, x_corr.shape

    if n < 0:
        mono = checkmonotony(x_corr)
        print 'mono after?',mono
    else:
        mono = checkmonotony(x_corr[:,sortidx])
        print 'mono after?',mono
    
    
    return x_corr
                
    
    
#    if n < 0:
        

    
    delta = 0.

    for i in np.arange(x.shape[0]-1):
        delta_tmp = x[i+1] - x[i]
	x[i] += delta
        
        if (delta_tmp < 0. and np.absolute(delta_tmp)>0.0003):
	    delta += np.absolute(delta_tmp)			
            
    x[-1] += delta
            
    return x
