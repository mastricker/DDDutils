#!/usr/bin/env python
"""
author: Markus Stricker

Functions written for the analysis of <eps_ij_000XXXXX.csv> files, which
come out of disloc when choosing option 19 and give the plastic strain
in voxels according to their swept area.
"""

import numpy as np

def eps_ij_csv_to_tensor(data, idxs):
    """
    Function for converting line based data from disloc option 19
    to tensor (np.array, 3x3).

    Args:
        data (np.array): loadtxt version of disloc output eps_ij_XXX.csv
        idxs (list(int)): contains number of slices in x/y/z dir [nx,ny,nz]

    Returns:
        data (np.array): same size, except column based eps_ij is now a 
            tensor.
        dx/dy/dz (int): extent of each slice
    """

    flag_idxs = False

    eps_tensor = np.empty([idxs[0],idxs[1],idxs[2],3,3])

    # Get dx/dy/dz for indexing if necessary
    if any (i > 1 for i in idxs):
        dx = data[0,0]*2.
        dy = data[0,1]*2.
        dz = data[0,2]*2.
        flag_idxs = True

    for i in range(data.shape[0]):
        if flag_idxs:
            idx = int(data[i,0]/dx)
            idy = int(data[i,1]/dy)
            idz = int(data[i,2]/dz)
        else:
            idx = idy = idy = 0

        eps = np.array([ [data[i,3],data[i,6],data[i,7]],
                         [data[i,6],data[i,4],data[i,8]],
                         [data[i,7],data[i,8],data[i,5]] ])

        eps_tensor[idx,idy,idz,:,:] = eps

    return eps_tensor, [dx,dy,dz]

def eps_ij_csv_to_tensor_single(data):
    """
    Function for converting line based data from disloc option 19
    to tensor (np.array, 3x3).

    Args:
        data (np.array): loadtxt version at stripped of all unnecessary values;
            size should be a line with 6 entries (xx, yy, zz, xy, xz, yz)

    Returns:
        data (np.array): same size, except column based eps_ij is now a 
            tensor.
    """

    eps_tensor = np.array([ [data[0], data[3], data[4]],
                            [data[3], data[1], data[5]],
                            [data[4], data[5], data[2]] ])

    return eps_tensor

def vMises_eq_strain(eps_ij):
    """
    Returns the von Mises average strain for the given strain tensor.

    Args:
        eps_ij: array 3x3
            Strain tensor with matching dimensions.

    Returns:
        eps_vM: float
            Von Mises equivalent strain.
    """
    
    eps_dev = eps_ij - 1./3.*np.trace(eps_ij)*np.identity(3)

    eps_vM = np.sqrt(2./3.*np.tensordot(eps_dev, eps_dev))

    return eps_vM
    
def eps_ij_to_eps_v_mises_eq(data):
    """
    Function calculates v Mises average strain for a given input size. 
    It expects the input to be the following:

    Args:
        data (np.array [ix,iy,iz,3,3]): Voxels should be indexed in three
            dimensions and eps_pl should be in tensor form.

    Returns:
        res (np.array [ix,iy,iz]): Von Mises equivalent plastic strain.
            one value per voxel
    """

    if not data.shape[-2:] == (3,3):
        print 'Can not process this type: data.shape', data.shape
        print 'data.shape[-2:]'.data.shape[-2:]
        exit

    N = len(data.shape)-2

    if N < 0 or N > 3:
        print 'Shape does not conform with function.'
        exit

    res = np.empty(data.shape[:-2])

    for i,_ in np.ndenumerate(res):
        eps = data[i]
        res[i] = vMises_eq_strain(eps)

    return res


def read_eps_pl_tensor(path):
    """
    Returns the available eps_ij_pl which has to be calculated before as a
    tensorial array.
    Each slip system has its own time series of the plastic strain tensor.

    The inputs, which are read here have to be prepared using

    ~/DDD/scripte/calculate_plastic_strain_macro_gs.py

    Args: 
        path (str): path to simulation folder

    Returns:
        eps_tensor (np.array): size is corresponing to [N,12,3,3], where N is
            the simulation steps, 12 is the number of slip systems and the 3x3
            rest ist the actual tensor

        time (np.array): corresponding time per data series
    """

    # load all time-series entries, splitted into contributions
    eps_xx = np.loadtxt(path+'/results/eps_ij_gs_macro_all/eps_gs_macro_xx.csv')
    eps_yy = np.loadtxt(path+'/results/eps_ij_gs_macro_all/eps_gs_macro_yy.csv')
    eps_zz = np.loadtxt(path+'/results/eps_ij_gs_macro_all/eps_gs_macro_zz.csv')
    eps_xy = np.loadtxt(path+'/results/eps_ij_gs_macro_all/eps_gs_macro_xy.csv')
    eps_xz = np.loadtxt(path+'/results/eps_ij_gs_macro_all/eps_gs_macro_xz.csv')
    eps_yz = np.loadtxt(path+'/results/eps_ij_gs_macro_all/eps_gs_macro_yz.csv')

    # and convert them to tensors
    N = eps_xx.shape[0] # number of entries

    time = eps_xx[:,1]

    eps_tensor = np.zeros([N,12,3,3])

    for i in range(N):
        for j in range(3,15):
            eps = np.array([ [eps_xx[i,j], eps_xy[i,j], eps_xz[i,j]],
                             [eps_xy[i,j], eps_yy[i,j], eps_yz[i,j]],
                             [eps_xz[i,j], eps_yz[i,j], eps_zz[i,j]]])
            eps_tensor[i,j-3,:,:] = eps

    return eps_tensor, time
