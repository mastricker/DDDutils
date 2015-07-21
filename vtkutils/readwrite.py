import vtk

"""
Helper functions to read and write vtk files
"""

# read vtkfile, return data + info for visual feedback
def read(filename):
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    header = reader.GetHeader()
    no_nodes = data.GetNumberOfPoints()
    no_connections = data.GetNumberOfCells()
    return (data, header, no_nodes, no_connections)

def write(data,filename):
    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(filename + ".vtk")
    writer.SetInput(data)
    writer.Write()
    print "Data written to file:\t " + filename +".vtk"
