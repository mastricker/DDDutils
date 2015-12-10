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

    eps_tensor = np.empty([idxs[0],idxs[1],idxs[2],3,3])

    # Get dx/dy/dz for indexing
    dx = data[0,0]*2.
    dy = data[0,1]*2.
    dz = data[0,2]*2.

    for i in range(data.shape[0]):
        idx = int(data[i,0]/dx)
        idy = int(data[i,1]/dy)
        idz = int(data[i,2]/dz)

        eps = np.array([ [data[i,3],data[i,6],data[i,7]],
                         [data[i,6],data[i,4],data[i,8]],
                         [data[i,7],data[i,8],data[i,5]] ])

        eps_tensor[idx,idy,idz,:,:] = eps

    return eps_tensor, [dx,dy,dz]
    
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

    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]

    res = np.empty([nx,ny,nz])
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):

                eps = data[i,j,k,:,:]

                eps_dev = eps - 1./3.*np.trace(eps)*np.identity(3)

                eps_eq = np.sqrt(2./3.*np.tensordot(eps_dev,eps_dev))

                res[i,j,k] = eps_eq

    return res
