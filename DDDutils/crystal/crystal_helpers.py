#!/usr/bin/env python
"""
author: Markus Stricker

Helper functions to deal with crystallographic rotations.
"""

import numpy as np

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
