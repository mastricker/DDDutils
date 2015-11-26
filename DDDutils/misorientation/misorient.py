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
import ..crystal 

# crystal system, fcc


def init_crystal_system(e1,e2,e3):
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

        m1 = float(m[0])
        m2 = float(m[1])
        m3 = float(m[2])

    # labels and directions to project misorientations on
    label0 = '_noproj'

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
