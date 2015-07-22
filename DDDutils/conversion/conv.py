import vtk
import networkx as nx
import numpy as np
import math

"""
Helper functions to change data from one format to another
"""

def nxtovtk(G):

    """
    BE PREPARED: THIS TAKES A VEEEEERY LONG TIME, EVEN FOR SMALL DATASETS
    input: networkx MultiGraph
    output: vtk data, POINTS + CELLS
    nodes attributes: position
    edge attributes: distance between nodes
    cell_data color scalars

    !!!-----------------------------------------------------------
    rearrange indices before converting to vtk!!!!!!!!!!!!!!!!!!!!
    !!!-----------------------------------------------------------
    """

    
    print "Converting graph to vtk-object:"
    # instanciate vtk object
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(G.nodes()))
    #print G.nodes()
    nodeids = np.array([[0,G.nodes()[0]]])
#    colors = vtk.vtkUnsignedCharArray()
#    colors.SetNumberOfComponents(4)
#    colors.SetName(seg_colors)

    # get new node ids for consecutive node numbering
    print "\t Storing new node Id's..."
    for i in range(len(G.nodes())):
        nodenr = G.nodes()[i]
        points.InsertPoint(i, (G.node[nodenr]['x'], G.node[nodenr]['y'], G.node[nodenr]['z']))
        if i == 0:
            pass
        else:
            new_col = np.array([i,G.nodes()[i]])
            nodeids = np.vstack((nodeids,new_col))
    #print nodeids
    grid.SetPoints(points)
    print "\t\t Done."
    # convert edge data to lines and add to grid
    # with checking for new node ids to keep the network data consistent
    print "\t Remapping edges with new Id's..."
    for m in G.edges():
        id_list=vtk.vtkIdList()
#        print m
        for nodeidold in m:
            n = 0
            while n <= range(len(nodeids[:,1])):
                if nodeidold == nodeids[n,1]:
                    nodeidnew = nodeids[n,0]
                    break
                else:
                    n += 1               
#            print nodeidnew
            id_list.InsertNextId(nodeidnew)
            cell_type = 3
        grid.InsertNextCell(cell_type, id_list)
        #print "added edge"
        id_list.Reset()
    print "\t\t Done."
    print "Conversion done."
    return grid

def nxtovtk_allpoint_only_edges(G_orig,G_clipped):
    """
    Write VTK file with all points, but only remaining edges.
    """
    
    print "Converting graph to vtk-object:"
    # instanciate vtk object
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(G_orig.nodes()))

    i = 0

    map_ids = np.empty([100000, 1])
    
    print "\t Storing new node Id's..."
    for n in G_orig.nodes():
        nodenr = n

        x = float(G_orig.node[nodenr]['x'])
        y = float(G_orig.node[nodenr]['y'])
        z = float(G_orig.node[nodenr]['z'])

        points.InsertPoint(i, (x, y, z))
        map_ids[n] = i

        
        i += 1

    # set all points from original dataset as new points
    grid.SetPoints(points)

    # Allocate storage for cells
    grid.Allocate(len(G_clipped.edges()))
    
    # now add edges from changed structure
    
    id_list=vtk.vtkIdList()
    id_list.Reset()

    for m in G_clipped.edges():
        for nodeid in m:
            mapped_id = map_ids[nodeid]
            id_list.InsertNextId(mapped_id)

        cell_type = 3
        grid.InsertNextCell(cell_type, id_list)
        id_list.Reset()

    print "\t\t Done."
    print "Conversion done. Number of cells:", grid.GetNumberOfCells()

    return grid


def vtktonx(data):
    data_nz = cpnonzerocolor(data)
    G = nx.Graph() # change to multigraph if structure needs to be consistent

    #nodeids = np.array([[0,data.GetCell(0).GetPointId(0)]])
    #print "Point id0 of cell 0:\t", data.GetCell(0).GetPointId(0)
    #idcntr = 1
    normal_colors = data_nz.GetCellData().GetScalars()
    ar0 = data_nz.GetCellData().GetArray(0) # coloring for normals
    ar1 = data_nz.GetCellData().GetArray(1) # dislocation character
    ar2 = data_nz.GetCellData().GetArray(2) # burgers index
    ar3 = data_nz.GetCellData().GetArray(3) # normal index
#    ar4 = data.GetCellData().GetArray(4) # cell ID
    
    for i in range(data_nz.GetNumberOfCells()):
        #print "colored:\t", i
#        print normal_colors.GetTuple(i)
#        print "0", ar0.GetTuple(i)
#        print "1", ar1.GetTuple(i)
#        print "2", ar2.GetTuple(i)
#        print "3", ar3.GetTuple(i)
        
        cell = data_nz.GetCell(i)

        id0 = cell.GetPointId(0)
        p1 = data_nz.GetPoint(id0)
        id1 = cell.GetPointId(1)
        p2 = data_nz.GetPoint(id1)
        distSquared = vtk.vtkMath.Distance2BetweenPoints(p1,p2)
        distance = math.sqrt(distSquared)
        G.add_node(int(id0))
        G.node[int(id0)]['x'] = float('%.4f' % p1[0])
        G.node[int(id0)]['y'] = float('%.4f' % p1[1])
        G.node[int(id0)]['z'] = float('%.4f' % p1[2])
        G.add_node(int(id1))
        G.node[int(id1)]['x'] = float('%.4f' % p2[0])
        G.node[int(id1)]['y'] = float('%.4f' % p2[1])
        G.node[int(id1)]['z'] = float('%.4f' % p2[2])
        G.add_edge(int(id0),int(id1), length = distance,
                   normal_color = ar0.GetTuple(i),
                   disloc_char = ar1.GetTuple(i),
                   burgers_index = ar2.GetTuple(i),
                   normal_index = ar3.GetTuple(i))
    #print G.nodes()
    print "...converted vtk to nx.Graph"
    #print nodeids
    return G


# copy cells with nonzero color to avoid junctions
def cpnonzerocolor(data):
    CellIdList = vtk.vtkIdList()
    CellIdList.Reset()
    colors = data.GetCellData().GetScalars()
    count = 0
    # Get indices from nonzero colored cells
    for i in range(data.GetNumberOfCells()):
        if ( colors.GetTuple(i) != (0.0, 0.0, 0.0, 0.0) ):
            CellIdList.InsertId(count, i)
            count = count + 1
#            print "Number Of cell:\t",i
#            print "Number of cell added:\t",count
        else:
            pass
    # Extract cells with nonzero color from data
    extractor = vtk.vtkExtractCells()
    extractor.SetInputData(data)
    extractor.SetCellList(CellIdList)
    extractor.Modified()
    extractor.Update()
    print "Number of nodes in clipped subvolume:\t",extractor.GetOutput().GetNumberOfPoints()
    print "Number of egdes in clipped subvolume:\t",extractor.GetOutput().GetNumberOfCells()
    # rearrange the id's of the cells and points - consecutive increasing ids

    print "...extracted zerocolor cells."
    
    return extractor.GetOutput()
