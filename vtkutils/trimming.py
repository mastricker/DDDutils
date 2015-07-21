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
