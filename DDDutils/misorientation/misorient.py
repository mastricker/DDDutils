#!/usr/bin/env python
"""
author: Markus Stricker

Function written for the analysis of <pos_deriv_00000050.dat> files, which
come out of disloc when choosing option 4.
Use a high discretization, when running disloc (min 20x20x20 voxels,
depending on sample size). Choose the same discretization in in x and z
direction - or change script. This script assumes a square cross-section 
with the same number of voxels in each direction.

Hint: EBSD measurements usually yield an average over approx. 100 nm.
But: It does not really make sense to choose a much larger discretization
than the average dislocation spacing. Try several discretizations for raw input


"""

import numpy as np
# relative imports
import DDDutils.crystal.crystal_helpers as ch

# crystal system, fcc


def init_crystal_system():
    """
    Return the orientations of the fcc lattice, according to the
    orientation of the specimen. These directions can be used to calculate
    misorientation along small angle grain boundaries.


    Args:
        e1,e2,e3 (double): Euler angles of crystallographic orientation
    
    Returns:
        labels (list of strings): labels for file-naming 
        orientations (list of np.arrays): x/y/z, burgersvectors and
            normals of the speciment, rotated by the euler angles
    """

    # Get sample orientation

    with open("../../normals_and_burgers.dat","r") as f:
        m = f.readlines()[0].split()[1:4]

    m1 = ch.grad_to_rad(float(m[0]))
    m2 = ch.grad_to_rad(float(m[1]))
    m3 = ch.grad_to_rad(float(m[2]))

    ME = ch.get_euler_matrix(m1,m2,m3)

    # labels and directions to project misorientations on

    label11 = '_100'
    label12 = '_010'
    label13 = '_001'

    label21 = '_110'
    label22 = '_101'
    label23 = '_011'
    label24 = '_-110'
    label25 = '_-101'
    label26 = '_0-11'

    label31 = '_1-1-1'
    label32 = '_111'
    label33 = '_-1-11'
    label34 = '_-11-1'

    labels = [label11, label12, label13,
              label21, label22, label23, label24, label25, label26,
              label31, label32, label33, label34]

    
    orientation11 = np.array([ 1., 0., 0.])
    orientation12 = np.array([ 0., 1., 0.])
    orientation13 = np.array([ 0., 0., 1.])

    orientation21 = np.array([-1., 0., 1.])
    orientation22 = np.array([ 0., 1., 1.])
    orientation23 = np.array([-1., 1., 0.])
    orientation24 = np.array([ 1., 1., 0.])
    orientation25 = np.array([ 0., 1.,-1.])
    orientation26 = np.array([ 1., 0., 1.])

    orientation31 = np.array([ 1.,-1.,-1.])
    orientation32 = np.array([ 1., 1., 1.])
    orientation33 = np.array([-1.,-1., 1.])
    orientation34 = np.array([-1., 1.,-1.])

    orientations = [orientation11, orientation12, orientation13,
                    orientation21, orientation22, orientation23,
                    orientation24, orientation25, orientation26,
                    orientation31, orientation32, orientation33,
                    orientation34]
    
    orientations = [vec/np.linalg.norm(vec) for vec in orientations]

    orientations_rot = [np.dot(ME,vec) for vec in orientations]

    return labels,orientations_rot

def index_raw_deriv_data(data):
    """
    Convert raw data from disloc pos_deriv* to indexed array.
    The choice which components of the derivatives is also made here.
    If you want to include the fem elasic part/ disloc elastic part/
    disloc plastic part

    Args:
        data (np.array): Raw data from DDD

    Returns:
        deriv_idx (np.array): indexed raw data
    """

    data_x = data[:,0]			
    data_y = data[:,1]
    data_z = data[:,2]
    data_x = data_x.tolist()		# convert numpy array to python list
    data_y = data_y.tolist()
    data_z = data_z.tolist()
    data_x = list(set(data_x))		# removing duplicates in List
    data_y = list(set(data_y))
    data_z = list(set(data_z))
    data_x = sorted(data_x)		# sort list
    data_y = sorted(data_y)
    data_z = sorted(data_z)
    
    idx1 = len(data_x)			
    idy1 = len(data_y)
    idz1 = len(data_z)

    deriv_idx = np.zeros([idx1,idy1,idz1,9])
    
    delta_x = data_x[1] - data_x[0]		# calculate step size
    delta_y = data_y[1] - data_y[0]             # for indexing
    delta_z = data_z[1] - data_z[0]

    for i in range(data.shape[0]):
        idx = int(round(  data[i,0] / delta_x)) # calculate current index
        idy = int(round( (data[i,1] - data_y[0]) / delta_y) )# correct offset of y-data
        idz = int(round(  data[i,2] / delta_z) )
        
        du_fem = np.zeros(9) # zeros for elastic part
#        du_fem = data[i,3:12] # comment in for inclusion of elastic part from fem

        du_del = np.zeros(9)
        du_del = data[i,12:21]
        du_dpl = np.zeros(9)
        du_dpl = data[i,21:30]

        du_tot = du_fem + du_del + du_dpl

        deriv_idx[idx,idy,idz,:] = du_tot

    return deriv_idx

def get_rotation_vector_voxel(data):
    """
    Calculate the rotation vector as a mean value for each voxel from 
    derivatives.

    Args:
        data (np.array): Indexed displacement derivatives for each voxel.

    Returns:
        w_voxel (np.array): Rotation vector for each voxel from average
        displacement gradients.
    """

    w_voxel = np.empty((data.shape[0]-1,data.shape[1]-1,data.shape[2]-1,4))

    for i in range(data.shape[0]-1):
        for j in range(data.shape[1]-1):
            for k in range(data.shape[2]-1):

                # get average displacement gradients per voxel
                du_sum = np.zeros([1,9])

                du_sum += data[i  ,j  ,k  ,:]
                du_sum += data[i  ,j+1,k  ,:]
                du_sum += data[i  ,j  ,k+1,:]
                du_sum += data[i  ,j+1,k+1,:]
                du_sum += data[i+1,j  ,k  ,:]
                du_sum += data[i+1,j+1,k  ,:]
                du_sum += data[i+1,j  ,k+1,:]
                du_sum += data[i+1,j+1,k+1,:]
                
                du_sum /= 8.

                w = np.empty(3)

                # crossproduct
                w[0] = du_sum[0,8] - du_sum[0,7]
                w[1] = du_sum[0,5] - du_sum[0,6]
                w[2] = du_sum[0,4] - du_sum[0,3]

                w = 0.5*w

                nu = np.linalg.norm(w)
                w_norm = w/nu

                w_voxel[i, j, k, :] = np.array([w_norm[0], w_norm[1], w_norm[2], nu])

    return w_voxel

def get_omega_voxel(data_w):
    """
    Calculates the rotation matrix based on the rotation vector.

    Args:
        data_w (np.array): Voxel based rotation vector.

    Returns:
        omega_voxel (np.array): Voxel indexed rotation matrices.
    """

    idx = data_w.shape[0]
    idy = data_w.shape[1]
    idz = data_w.shape[2]
  
    omega_voxel = np.empty([idx, idy, idz, 3, 3])
  
    for i in range(idx):
        for j in range(idy):
            for k in range(idz):

	        Omega = np.ones((3,3))
                r = data_w[i,j,k,0:3]
        
	        nu = data_w[i,j,k,3]

                rot_matrix = np.eye(3)*np.cos(nu) \
                             + (1-np.cos(nu))*np.outer(r,r) \
                             + np.sin(nu) \
                             * np.array ([[0., r[2], -r[1]],[-r[2], 0, r[0]],[r[1], -r[0], 0.]])

                omega_voxel[i,j,k,:,:] = rot_matrix

    return omega_voxel

def calculate_gamma(R1, R2):
    """
    Calculates the misorientation angle between two given rotation
    matrices.

    Args:
        R1,R2 (np.array): 3x3 rotation matrix

    Returns:
        g (double): misorientation angle in rad
    """

    Rmis = np.dot(R2,np.linalg.inv(R1))

    val = ( np.trace(Rmis) - 1. ) / 2.

    # val < -1 ----> -1
    if (val > 1.0):
        print 'case1',val
        val = 1.0
    elif( val < -1.0):
        print 'case2',val
        val = -1.0
        
    g = np.arccos(val)

    if (np.isnan(g)):
        g = 0.
        print 'should not occur: nan error in calculate_gamma'

    return g
    
def misorientation_wrt_coord(omega,coords):
    """
    Returns the misorientation of each voxel with respect to a
    reference coordinate.

    Args:
        omega (np.array): Voxel based rotation matrices
        coords (list(int)): 3 indices in a list

    Returns:
        misorient (np.array): Misorientation angle of each voxel
            with respect to given reference voxel.
    """

    idx = Omega.shape[0]
    idy = Omega.shape[1]
    idz = Omega.shape[2]
  
    misorient = np.zeros([idx, idy, idz])
  
    Rref = omega[coords[0],coords[1],coords[2],:,:]

    for i in range(idx):
        for j in range(idy):
            for k in range(idz):
	        R = Omega[i,j,k,:,:]
	
	        gamma, flag_nan = calculate_gamma(Rref, R)
	
	        if not flag_nan:
	            misorient[i,j,k] = rad_to_grad(gamma)
	        else:
                    print 'set gamma = 0'
	            misorient[i,j,k] = 0.
	  
  return misorient
