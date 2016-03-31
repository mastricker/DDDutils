#!/usr/bin/env python
"""
author: Markus Stricker

Helper functions to deal with crystallographic operations.
"""

import numpy as np


# helper, needed for some of the functions below
n2b = [  [1,4],  [1,5],  [1,6],\
         [2,1],  [2,3],  [2,5],\
         [3,2],  [3,3],  [3,6],\
         [4,1],  [4,2],  [4,4]  ]


map_schmidboas = [ 2, 3, 1, 6, 4, 5, 7, 9, 8, 11, 10, 12]

gsystems = ("A2", "A3", "A6",\
            "B2", "B4", "B5",\
            "C1", "C3", "C5",\
            "D1", "D4", "D6")


################### FUNCTIONS

def grad_to_rad(alpha_deg):
    """ Convert degrees to radians. """
    
    alpha_rad = (alpha_deg/180.)*np.pi
    return alpha_rad

def rad_to_grad(alpha_rad):
    """ Convert radians to degrees. """
    
    alpha_deg = (alpha_rad*180.)/np.pi
    return alpha_deg

def get_euler_matrix(phi_rad, theta_rad, psi_rad):
    """
    Returns the Euler matrix according to the disloc convention

    Args:
        phi_rad, theta_rad, psi_rad (double): Euler angles in rad

    Returns:
        ME (np.array(), 3x3): Euler rotation matrix
    """
    
    D = np.array(([ np.cos(phi_rad), np.sin(phi_rad), 0.],
                  [-np.sin(phi_rad), np.cos(phi_rad), 0.],
                  [ 0., 0., 1.]))
    C = np.array(([ 1., 0., 0.],
                  [ 0., np.cos(theta_rad), np.sin(theta_rad)],
                  [ 0.,-np.sin(theta_rad), np.cos(theta_rad)]))
    B = np.array(([ np.cos(psi_rad), np.sin(psi_rad), 0.],
                  [-np.sin(psi_rad), np.cos(psi_rad), 0.],
                  [ 0., 0., 1.]))

    ME = np.dot(B,(np.dot(C,D)))

    return ME

def get_normals_and_burgers(path):
    """
    Returns the normals and burgers vectors of the given path.

    Args: 
        path (str): path to simulation folder

    Returns:
        normals, burgers (dict): Dictionary with n1/b2/etc as numpy arrays
    """


    normals = {}
    burgers = {}
    
    with open(path+'/normals_and_burgers.dat') as f:
        n = 0
        while(n<11):
            line = f.readline().split()
            if (n>0):
#            print('n={}, line: {}'.format(n,line))
                if (n<5):
                    normals['n'+str(line[0])] = np.array([float(line[1]), float(line[2]), float(line[3])])
                else:
                    burgers['b'+str(line[0])] = np.array([float(line[1]), float(line[2]), float(line[3])])
            n += 1

    return normals,burgers

def get_schmid_tensor(normals, burgers):
    """
    Function, which calculates the Schmid tensor:

    M = d \otimes n

    where d is the slip direction b_vec/b and n is the glide plane normal.

    Args:
        normals (dict): a dictionary contaning all normals, with names like n1, n2, n3
        
        burgers (dict): a dictionary containing all burgers vector, naming b1, b2, ...

    Returns:
        M_alpha (np.array): size is [12,3,3], 12 slip system, 3x3: Schmid tensor
    """

    M_alpha = np.zeros([12,3,3])

    for idx, (i,j) in enumerate(n2b):

        n = normals['n'+str(i)].copy()
        n = n/np.linalg.norm(n).copy()
        b = burgers['b'+str(j)]
        b = b/np.linalg.norm(b)

        M_alpha[idx,:,:] = np.outer(b,n)

    return M_alpha


def get_glissile_reaction_matrix():
    # glissile reaction matrix
    reaction_matrix = np.array([[( 3, 6),( 2, 9)],\
                                [( 1, 9),( 3,12)],\
                                [( 1, 6),( 2,12)],\
                                [( 6, 8),( 5, 2)],\
                                [( 6,10),( 4, 2)],\
                                [(10, 5),( 4, 8)],\
                                [( 9, 5),( 3, 8)],\
                                [( 3, 7),( 9,11)],\
                                [( 5, 7),( 8,11)],\
                                [( 7,12),( 1,11)],\
                                [( 1,10),( 4,12)],\
                                [( 7,10),( 4,11)] ])
    # shift index to python notation
    for i in reaction_matrix:
        for j in i:
            j[0] -= 1
            j[1] -= 1

    return reaction_matrix


def get_schmidboas_dict():
    # make dictionary
    dict_schmidboas = {}

    for i in range(12):
        dict_schmidboas[map_schmidboas[i]-1] = gsystems[i]

def get_schmid_factor(n,b,taxis=np.array([0.,1.,0.])):
    """
    Function for calculating the Schmid factor.

    Calculated by using cos(gamma)*cos(kappa).

    Where gamma is the angle between tensile axis and glide 
    direction and kappa the angle between tensile axis and glide
    plane normal.

    Args:
        n (np.array): normal of glide plane
    
        b (np.array): burgers vektor

        taxis (np.array): tensile axis (optional, assuming standard
            DDD tensile direction)

    Return:
        m (float): Schmid factor
    """

    b /= np.linalg.norm(b)
    n  /= np.linalg.norm(n)
    taxis /= np.linalg.norm(taxis)

    return np.dot(taxis,b)*np.dot(taxis,n)
    
