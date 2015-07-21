import vtk
import math
import numpy as np

"""
Helper functions for clipping vtk geometries
"""

def clip(data, extent):
#   print "CLIP CLIP CLIP CLIP CLIP CLIP CLIP CLIP ", extent
    boxclip = vtk.vtkBoxClipDataSet()
    boxclip.GenerateClipScalarsOn()
    #boxclip.GenerateClippedOutputOn()
    boxclip.SetBoxClip(*extent)
    boxclip.SetInputData(data)
    boxclip.Update()
    return boxclip.GetOutput()

# cut the volume with a line, return points in plane
# even if they are not specified in the original data
def cutplanes(data, origin, normal):
    # create cutting plane
    plane = createplane(origin,normal)
    # create cutter
    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(data)
    cutter.Update()
    #print cutter
    cutslice = cutter.GetOutput()
    npoints = cutslice.GetNumberOfPoints()
    ncells = cutslice.GetNumberOfCells()
    #print cutslice
    return cutslice, npoints, ncells
    

# create a plane to cut, here it cuts in the XZ plane
# yz normal=(1,0,0); xz = (0,1,0); xy =(0,0,1);
def createplane(origin, normal):
    plane = vtk.vtkPlane()
    plane.SetOrigin(origin) #(x,y,z)
    plane.SetNormal(normal) #(x,y,z)
    return plane

# create a cube
def createcube(box_extent):
    cube = vtk.vtkCubeSource()
    cube.SetBounds([0., box_extent[0], 0., box_extent[1], 0., box_extent[2]])
    cube.Update()

    return cube

# clip a volume with a normal and two points
def clipwithnormalandorigins(data,box,normal, p, delta):
    """
    http://www.vtk.org/doc/nightly/html/classvtkBoxClipDataSet.html
    SetBoxClip:( 	const double *  	n0,
		const double *  	o0,
		const double *  	n1,
		const double *  	o1,
		const double *  	n2,
		const double *  	o2,
		const double *  	n3,
		const double *  	o3,
		const double *  	n4,
		const double *  	o4,
		const double *  	n5,
		const double *  	o5 
	) 	

    assuming clipping orthogonal to y-axis
    """

    delta2 = delta/2.

    # x-direction
    p0 = p - delta2 * normal
    n0 = -normal
    p1 = p + delta2 * normal
    n1 = normal

    # y-direction
    p2 = np.array([p[0], 0., p[2]])
    n2 = np.array([0., -1., 0.])
    p3 = np.array([p[0], box[1], p[1]])
    n3 = -n2

    # z-direction
    p4 = np.array([box[0], p[1] ,p[2]])
    n4 = np.array([1., 0., 0.])
    p5 = np.array([0. ,p[1] ,box[2]])
    n5 = -n4
    
    extent = np.array([n0, p0, n1, p1, n2, p2, n3, p3, n4, p4, n5, p5])

    data_clip = clip(data,extent)

    return data_clip
    
