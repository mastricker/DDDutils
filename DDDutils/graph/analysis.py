import networkx as nx
import numpy as np
import math

"""
Functions to analyse graph structures from dislocation data
"""

######################################################################
# this function deletes all subgraphs from a set which
# do not reach from adjacent sites of volume boundaries
# defined by VOI of boxclip
######################################################################


def rmsubg(G, direction, extent):
    if direction == 'x':
        boundaries = [extent[0],extent[1]]
    elif direction == 'y':
        boundaries = [extent[2],extent[3]]
    elif direction == 'z':
        boundaries = [extent[4],extent[5]] 
    else:
        print "Direction not available."
        return
    print "Boundaries in ",direction,"-direction:\t", boundaries
    # get subgraphs
    components = nx.connected_component_subgraphs(G)
    
    # evaluate if subgraph has nodes in both planes of interest
    for i in range(len(components)):
        lowercounter = 0
        uppercounter = 0
        for n in components[i]:
#            print "#########"
#            print components[i].node[n][direction]
#            print direction, "-component of node ",n," in component ",i, " should be ", boundaries
#            print "and is ", components[i].node[n][direction]
#            print "------------------------------------------------------------------------------"
#            print "#########"
            if components[i].node[n][direction] == boundaries[0]:
                lowercounter += 1
            elif components[i].node[n][direction] == boundaries[1]:
                uppercounter += 1
            else:
#                print "node not @ boundary"
                pass
        if lowercounter >=1 and uppercounter >= 1:
            pass
        else:
            #print components[i].edges()
            G.remove_edges_from(components[i].edges())
            G.remove_nodes_from(components[i].nodes())
            
#        print "------------lower ctr:\t", lowercounter
#        print "------------upper ctr:\t", uppercounter
        #comp2  = nx.connected_component_subgraphs(G)
        #print"................................................................."
        #print len(comp2)
    print "...removed ", len(components)-len(nx.connected_component_subgraphs(G)), "/",len(components), "subgraphs"
    return G



######################################################################
# function for analysis of subgraphs I/O - flow
# inout is a list of in/out (lower/upper) connections in VOI
# return is ordered for subgraphs:
# inout = ([[subg1_in, subg1_out],
#           [subg2_in, subg2_out],
#          ..................... ])
# new inout_sum = ([[ins, outs]])
######################################################################

def iosubg(G, direction, extent):
    if direction == 'x':
        boundaries = [extent[0],extent[1]]
    elif direction == 'y':
        boundaries = [extent[2],extent[3]]
    elif direction == 'z':
        boundaries = [extent[4],extent[5]] 
    else:
        print "iosubg: Direction not available."
        return
    components = nx.connected_component_subgraphs(G)
    inout = np.zeros((len(components),2))
#    print "inout before"
#    print inout
    #delta = np.zeros((len(components),1))
    for i in range(len(components)):
        incntr = 0
        outcntr = 0
        for n in components[i]:
            if components[i].node[n][direction] == boundaries[0]:
                outcntr += 1
            elif components[i].node[n][direction] == boundaries[1]:
                incntr += 1
            else:
                pass
        inout[i,:] = [incntr, outcntr]
#    print "inout after"
#    print inout
    inout_sum = np.array([[inout[:,0].sum(), inout[:,1].sum()]])
#    print "iosubg: inout_sum:"
#    print inout_sum
    return inout_sum


        
######################################################################
# function for analysis of subgraphs mean distance between bounding box
# u: list of nodeids on lower bound
# v: list of nodeid on upper bound
######################################################################

def meanpathlengthsubg(G, direction, extent):
    if direction == 'x':
        boundaries = [extent[0],extent[1]]
    elif direction == 'y':
        boundaries = [extent[2],extent[3]]
    elif direction == 'z':
        boundaries = [extent[4],extent[5]] 
    else:
        print "meanpathlengthsubg: Direction not available."
        return
    if boundaries[1] == boundaries[0]:
        mean_pathlength_per_subg =  [0]
        no_of_paths_per_subg = [0]
        return mean_pathlength_per_subg, no_of_paths_per_subg
    else:
        components = nx.connected_component_subgraphs(G)
        mean_pathlength_per_subg=[]
        mean_pl_total=[]
        no_of_paths_per_subg=[]
        for i in range(len(components)):
            # build a list of nodes for path length calculation
            u = [] # list of nodes on lower bound
            v = [] # list of nodes on upper bound
            for n in components[i]:
                if components[i].node[n][direction] == boundaries[0]:
                    u.append(n)
                elif components[i].node[n][direction] == boundaries[1]:
                    v.append(n)
                else:
                    pass
            mean = []
            for unode in u:
                for vnode in v:
                    pathlength = 0
                    if nx.has_path(components[i], unode, vnode) == True:
                        path = nx.shortest_path(components[i],unode,vnode,'length') # optional 'length'                        
                        for node in range(len(path)-1):
                            """ 
                            attention: components is a dictionary of dictionaries
                            therefore extra [0] at the end before 'length'.    
                            """
                            pathlength += components[i][path[node]][path[node+1]]['length']
                        #print "recent pathlength: ------------->", pathlength
                        mean.append(pathlength)
                    else:
                        pass
            print "#no of shortest paths in subg ",i,":\t", len(mean)
            for pl in mean:
                mean_pl_total.append(pl)
            no_of_paths_per_subg.append(len(mean))
            mean_pathlength_per_subg.append(sum(mean)/len(mean))
    return mean_pathlength_per_subg, no_of_paths_per_subg, mean_pl_total



######################################################################
# function for analysis of subgraphs mean distance between bounding box
# u: list of nodeids on lower bound
# v: list of nodeid on upper bound
######################################################################

def pathlength(G, direction, extent):
    """
    Parameters
    ----------
    G : standard graph
        with vertices (attributes are coordinates) and edges
        (attributes are length)
    
    direction : string, {x, y, z}
        for determining in which direction the shortest path search
        shall be performed
        
    extent : np.array [1x6] - xmin, xmax, ymin, ymax, zmin, zmax
        coordinates for a box which gives the extent of VOI
        from partitioning in main function
        
    Returns
    ------- 
    mean_pl : mean pathlength, float
        normalized to delta-direction
        over all paths through VOI
    
    std_dev : standard deviation, float
        for all paths through VOI
        
    no_paths : int
        number of paths in subsample in given direction through VOI
    """
    
    # get boundaries from direction
    if direction == 'x':
        boundaries = [extent[0],extent[1]]
    elif direction == 'y':
        boundaries = [extent[2],extent[3]]
    elif direction == 'z':
        boundaries = [extent[4],extent[5]] 
    else:
        print "pathlength: Direction not available."
        return
    
    # if box is just a plane, return this:
    if boundaries[1] == boundaries[0]:
        mean_pl =  [0]
        std_dev = [0]
        no_pls = [0]
        return mean_pl, std_dev, no_pls
    
    # for non-plane extents:
    else:
        # get delta for normalizing pathlengths to subsample extent
        delta = math.fabs(boundaries[1]-boundaries[0])
        
        # get independent subgraphs
        components = nx.connected_component_subgraphs(G)
        
        # initialize return values
        pls = []
        no_pls=[]
        
        # go through independent subgraphs and look for shortest paths
        for i in range(len(components)):
            """
            build a list of nodes for path length calculation
            which lie either on lower or upper bound of VOI
            """
            u = [] # list of nodes on lower bound
            v = [] # list of nodes on upper bound
            for n in components[i]:
                if components[i].node[n][direction] == boundaries[0]:
                    u.append(n)
                elif components[i].node[n][direction] == boundaries[1]:
                    v.append(n)
                else:
                    pass
            # collector for shortest paths found
            pls_coll = []
            for unode in u:
                for vnode in v:
                    pathlength = 0
                    if nx.has_path(components[i], unode, vnode) == True:
                        path = nx.shortest_path(components[i],unode,vnode,'length') # optional 'length'                        
                        for node in range(len(path)-1):
                            """
                            ! obsolete with standard graph:
                            attention: components is a dictionary of dictionaries
                            therefore extra [0] at the end before 'length'.    
                            """
                            pathlength += components[i][path[node]][path[node+1]]['length']
                        #print "recent pathlength: ------------->", pathlength
                        # normalize new found shortest path and add to collection
                        pls_coll.append(pathlength/delta)
                    else:
                        pass
            print "#no of shortest paths in subg ",i,":\t", len(pls_coll)
            for pl in pls_coll:
                pls.append(pl)
            no_pls.append(len(pls_coll))
        # calculate mean of all paths
        pls = np.array([pls])
        mean_pl = pls.mean()
        std_dev = np.std(pls)
        no_pls  = sum(no_pls)
    return mean_pl, std_dev, no_pls




######################################################################
# function for analysis of subgraphs: hotspots, degree
# clustering coefficient (of a node): is the fraction of pairs of its neighbors
#                                     that have edges between one another
# hub: maximum degree of a node in subgraph
######################################################################

def find_all_paths(graph, start, end):
    path  = []
    paths = []
    queue = [(start, end, path)]
    while queue:
        start, end, path = queue.pop()
        print 'PATH', path

        path = path + [start]
        if start == end:
            paths.append(path)
        for node in set(graph[start]).difference(path):
            queue.append((node, end, path))
    return paths
                
            
        
######################################################################
# function for analysis of subgraphs: hotspots, degree
# clustering coefficient (of a node): is the fraction of pairs of its neighbors
#                                     that have edges between one another
# hub: maximum degree of a node in subgraph
######################################################################

def findhotspotssubg(G, direction, extent):
    if direction == 'x':
        boundaries = [extent[0],extent[1]]
    elif direction == 'y':
        boundaries = [extent[2],extent[3]]
    elif direction == 'z':
        boundaries = [extent[4],extent[5]] 
    else:
        print "findhotspotssubg: Direction not available."
        return
    if boundaries[1] == boundaries[0]:
        pass
    else:
        components = nx.connected_component_subgraphs(G)
        hubs = []
        clust = []
        for i in range(len(components)):
            # convert subgraph from MultiGraph to Graph 
            nodes = components[i].nodes()
            edges = components[i].edges()
            h = nx.Graph()
            h.add_nodes_from(nodes)
            h.add_edges_from(edges)
            # get degree and clustering
            degree_sequence = nx.degree(h).values()
            hub = max(degree_sequence)
            hclust = nx.clustering(h).values()
            hclust_max = max(hclust)
            hubs.append(hub)
            clust.append(hclust_max)
        return hubs, clust
    


######################################################################
# function for analysis of subgraphs: fraction of volume in dislocation
# pipes
# see Hart [1957]: On the Role of Dislocations in Bulk Diffusion
# looking for 'f'
######################################################################

def pipevolfraction(G, lattice_constant, a, extent):
        """
        Parameters:
        -----------
        G : nx.Graph()
            Graph, in this case imported from vtk. Edge attributes have
            to be the distance between the nodes! Look into vtk2graph.py
            for how the vtk cells are converted into graph edges.
            
        lattice_constant : 6*float
            This number is the lattice constant to scale the atomic distances
            in the graph edge information and the 
            
        a : float
            This value is defined in main.py and represents the effective
            dislocation pipe radius in [m] 
            
        extent : 6*float
            This vector, in lattice units gives the surrounding box of the
            given graph. It is used to calculate the fraction of volume.
            
        Return:
        -------
        f : float
            Fraction of volume which is in pipes. Unit: [m]
        """
        
        pipe_lengths = nx.get_edge_attributes(G, 'length')
        
        l_pipes = 0
        
        for i,j in pipe_lengths:
            l_pipes += pipe_lengths[(i,j)]
        
        # scale length of pipe network to [m]
        l_pipes = lattice_constant * l_pipes
        
        # get volume of pipes: length*cross_section
        v_pipes = l_pipes*np.pi*a**2
        
        # get complete volume of extent, in [m]
        delta_x = lattice_constant*(extent[1] - extent[0])
        delta_y = lattice_constant*(extent[3] - extent[2])
        delta_z = lattice_constant*(extent[5] - extent[4])
        
        vol = delta_x * delta_y * delta_z
        
        # fraction of volume in dislocation pipes
        f = v_pipes / vol
        
        return f
