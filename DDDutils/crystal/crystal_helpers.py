#!/usr/bin/env python
"""
author: Markus Stricker

Helper functions to deal with crystallographic operations.
"""

import numpy as np
import math

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

gp_normal_ref = np.array([ [ 1., -1., -1.],\
                           [ 1.,  1.,  1.],\
                           [-1., -1.,  1.],\
                           [-1.,  1., -1.] ])

burger_ref = np.array([ [-1., 0., 1.],\
                        [ 0., 1., 1.],\
                        [-1., 1., 0.],\
                        [ 1., 1., 0.],\
                        [ 0., 1.,-1.],\
                        [ 1., 0., 1.] ])

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


def get_euler_angles(pathtograin):
    """
    Function for reading the Euler angles of a simulation.

    Args: 
        pathtograin (str): path to simulation folder and desired grain input file.

    Returns:
        e1,e2,e3 (float): Euler angles from path.
    """

    e1 = e2 = e3 = 0.
    with open(pathtograin, 'r') as fin:
        for i, line in enumerate(fin):
            if i==0:
                e1 = float(line.strip().split()[1])
                e2 = float(line.strip().split()[2])
                e3 = float(line.strip().split()[3])


    return e1, e2, e3

    
    
def calculate_rotation_matrix_lab_sys(proj_x, proj_z):
    """
    Function for calculating the rotation matrix. Assumption:
    New axis for x and z are given.
    """

    e1 = np.array([1.,0.,0.])
    e2 = np.array([0.,1.,0.])
    e3 = np.array([0.,0.,1.])

    g1 = proj_x / np.linalg.norm(proj_x)
    g3 = proj_z / np.linalg.norm(proj_z)

    g2 = np.cross(g3,g1)
    g2 /= np.linalg.norm(g2)

    MEuler = np.array(([np.dot(g1,e1),np.dot(g1,e2),np.dot(g1,e3)],
                       [np.dot(g2,e1),np.dot(g2,e2),np.dot(g2,e3)],
                       [np.dot(g3,e1),np.dot(g3,e2),np.dot(g3,e3)]))


    return MEuler

# def calculate_rotation_matrix_arbitrary(theta, u):
#     """
#     Function for calculating the rotation matrix of a an angle theta
#     around an arbitrary axis u.

#     Inspiration from:
#     http://stackoverflow.com/questions/17763655/
#     rotation-of-a-point-in-3d-about-an-arbitrary-axis-using-python

#     Args:
#         theta (float): angle in degrees

#         u (np.array, 3): rotation axis

#     Return:
#         Rotation matrix
#     """

#     theta = grad_to_rad(theta)
#     u /= np.linalg.norm(u)

#     return [[np.cos(theta) + u[0]**2 * (1-np.cos(theta)), 
#              u[0] * u[1] * (1-np.cos(theta)) - u[2] * np.sin(theta), 
#              u[0] * u[2] * (1 - np.cos(theta)) + u[1] * np.sin(theta)],
#             [u[0] * u[1] * (1-np.cos(theta)) - u[2] * np.sin(theta),
#              np.cos(theta) + u[1]**2 * (1-np.cos(theta)),
#              u[1] * u[2] * (1 - np.cos(theta)) + u[0] * np.sin(theta)],
#             [u[0] * u[2] * (1-np.cos(theta)) - u[1] * np.sin(theta),
#              u[1] * u[2] * (1-np.cos(theta)) - u[0] * np.sin(theta),
#              np.cos(theta) + u[2]**2 * (1-np.cos(theta))]]

def calculate_rotation_matrix_arbitrary(theta, axis):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    Inspiration:
    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    Args:
        theta (float): angle in degrees

        axis (np.array, 3): rotation axis

    Return:
        Rotation matrix
    """
    axis = np.asarray(axis)
    theta = np.asarray(grad_to_rad(theta))
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d

    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def get_intersection_lines(vec,gp_normals):
    
    vec_intersect = np.empty(gp_normals.shape)

    for i in range(gp_normals.shape[0]):
        gp_n = gp_normals[i,:]
        
        dir_intersect = np.cross(vec,gp_n)
        
        vec_intersect[i,:] = dir_intersect


    return vec_intersect


def get_rotation_matrix_with_proj(h1,g1,h2,g2):
    """
    h1: (new) reference x
    g1: direction to be projected on x

    h2: (new) reference y
    g2: direction to be projected on y
    """

    e1 = np.array([ 1., 0., 0.])
    e2 = np.array([ 0., 1., 0.])
    e3 = np.array([ 0., 0., 1.])

    h1 /= np.linalg.norm(h1)
    g1 /= np.linalg.norm(g1)

    h2 /= np.linalg.norm(h2)
    g2 /= np.linalg.norm(g2)

    h3 = np.cross(h1,h2)
    g3 = np.cross(g1,g2)

    MEuler_g_to_e = np.array(([np.dot(g1,e1),np.dot(g1,e2),np.dot(g1,e3)],
                          [np.dot(g2,e1),np.dot(g2,e2),np.dot(g2,e3)],
                          [np.dot(g3,e1),np.dot(g3,e2),np.dot(g3,e3)]))

    MEuler_h_to_e = np.array(([np.dot(h1,e1),np.dot(h1,e2),np.dot(h1,e3)],
                              [np.dot(h2,e1),np.dot(h2,e2),np.dot(h2,e3)],
                              [np.dot(h3,e1),np.dot(h3,e2),np.dot(h3,e3)]))
    
    MEuler_e_to_h = MEuler_h_to_e.transpose()

    h1_new0 = np.dot(MEuler_g_to_e, e1)
    h2_new0 = np.dot(MEuler_g_to_e, e2)
    h3_new0 = np.dot(MEuler_g_to_e, e3)
    
    h1_new = np.dot(MEuler_e_to_h, h1_new0)
    h2_new = np.dot(MEuler_e_to_h, h2_new0)
    h3_new = np.dot(MEuler_e_to_h, h3_new0)

    MEuler_h_new_to_e = np.array(([np.dot(h1_new,e1),
                                   np.dot(h1_new,e2),
                                   np.dot(h1_new,e3)],
                                  [np.dot(h2_new,e1),
                                   np.dot(h2_new,e2),
                                   np.dot(h2_new,e3)],
                                  [np.dot(h3_new,e1),
                                   np.dot(h3_new,e2),
                                   np.dot(h3_new,e3)]))

    M = MEuler_h_new_to_e.transpose()

    return M
        
