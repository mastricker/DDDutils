import vtk
from vtk.util import numpy_support as VN
import numpy as np

def convert_vtk_to_np_stress(data):
    """
    Function for changing point data of vtk file to stress tensor.
    
    Args:
        data (vtk): reader.GetOutput()

    Returns:
        sig (np.array): stress tensor per point 
    """

    n_nodes = data.GetNumberOfPoints()

    s_xx = VN.vtk_to_numpy(data.GetPointData().GetArray('s_tot_x'))
    s_yy = VN.vtk_to_numpy(data.GetPointData().GetArray('s_tot_y'))
    s_zz = VN.vtk_to_numpy(data.GetPointData().GetArray('s_tot_z'))
    s_xy = VN.vtk_to_numpy(data.GetPointData().GetArray('s_tot_xy'))
    s_xz = VN.vtk_to_numpy(data.GetPointData().GetArray('s_tot_xz'))
    s_yz = VN.vtk_to_numpy(data.GetPointData().GetArray('s_tot_yz'))
    

    sig_tensor = np.empty([n_nodes,3,3])

    for i in range(n_nodes):
        sig_tensor[i,0,0] = s_xx[i]
        sig_tensor[i,1,1] = s_yy[i]
        sig_tensor[i,2,2] = s_zz[i]
        sig_tensor[i,0,1] = s_xy[i]
        sig_tensor[i,1,0] = s_xy[i]
        sig_tensor[i,0,2] = s_xz[i]
        sig_tensor[i,2,0] = s_xz[i]
        sig_tensor[i,1,2] = s_yz[i]
        sig_tensor[i,2,1] = s_yz[i]

        
    return sig_tensor

