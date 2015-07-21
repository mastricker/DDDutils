import vtk

"""
Helper functions for clipping vtk geometries
"""

def clip(data, extent):
    print "CLIP CLIP CLIP CLIP CLIP CLIP CLIP CLIP ", extent
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
    cutter.SetInput(data)
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
