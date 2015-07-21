import vtk
import networkx as nx
import numpy as np

"""
Helper functions to change data from one format to another
"""

def nxtovtk(G):

    # convert graph labels to consecutive integers for vtk numbering
    # does not keep up with the connection data
    G = nx.convert_node_labels_to_integers(G, first_label=0)
    
    # instanciate vtk object
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(G.nodes()))
    
#    colors = vtk.vtkUnsignedCharArray()
#    colors.SetNumberOfComponents(4)
#    colors.SetName(seg_colors)
        
    # convert node data to points
    for n in G.nodes():
        points.InsertPoint(n, (G.node[n]['x'], G.node[n]['y'], G.node[n]['z']))
        #print n, G.node[n]['x'], G.node[n]['y'], G.node[n]['z']
    grid.SetPoints(points)

    # convert edge data to lines and add to grid
    for m in G.edges():
        id_list=vtk.vtkIdList()
#        print m
        for nodeid in m:
#            print nodeid
            id_list.InsertId(6, nodeid)
            cell_type = 3
        grid.InsertNextCell(cell_type, id_list)
        id_list.Reset()
    print "...converted nx.Graph to vtk object"
    return grid
