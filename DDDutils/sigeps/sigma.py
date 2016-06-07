"""
author: Markus Stricker

Functions written for analysis of stress data.
"""

import numpy as np


def von_Mises_stress(sigma):
    """
    Convert stress tensor to von Mises stress.

    Args:
        sigma (np.array [N,3,3]): stress data in tensor form with N entries

    Returns:
        sigma_vm (np.array [N,1]): von Mises stress per data point
    """

    n_points = sigma.shape[0]

    sigma_vm = np.empty([n_points])

    for i in range(n_points):
        sig_tmp = sigma[i,:,:]

        # sigma_vm[i] = np.sqrt(sig_tmp[0,0]**2 + sig_tmp[1,1]**2\
        #                    + sig_tmp[2,2]**2\
        #                    - sig_tmp[0,0] * sig_tmp[1,1]\
        #                    - sig_tmp[0,0] * sig_tmp[2,2]\
        #                    - sig_tmp[1,1] * sig_tmp[2,2]\
        #                    + 3 * (sig_tmp[0,1]**2 + sig_tmp[0,2]**2
        #                           + sig_tmp[1,2]**2 ))

        sigma_vm[i] = np.sqrt( 0.5 
                               * ( (sig_tmp[0,0] - sig_tmp[1,1])**2 
                                  +(sig_tmp[1,1] - sig_tmp[2,2])**2
                                  +(sig_tmp[2,2] - sig_tmp[0,0])**2
                                   + 6*(sig_tmp[0,1]**2 
                                      + sig_tmp[1,2]**2
                                      + sig_tmp[0,2]**2)))


    return sigma_vm
        
