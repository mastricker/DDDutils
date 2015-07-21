import numpy as np

# partitioning: BOTTOM TO TOP (btt)
# input:
# - sample_box = (xmin, xmax, ymin, ymax, zmin, zmax)
# - direction = 'x'/'y'/'z'
# - partition = int

def getpartitioningbtt(box, direction, partition, cutoff=0):
    if direction == 'x':
        cut = cutoff*(box[1]-box[0])
        delta = []
        upper = box[1]-cut
        lower = box[0]+cut
        deltax = (upper-lower)/partition
        partition_coordinates = np.array([[lower,lower,box[2],box[3],box[4],box[5]]])
        delta.append(deltax)
        for i in range(partition):
            new_row = np.array([lower,lower+(i+1)*deltax,box[2],box[3],box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append((i+1)*deltax)
        return partition_coordinates, delta
    elif direction == 'y':
        cut = cutoff*(box[3]-box[2])
        delta = []
        upper = box[3] - cut
        lower = box[2] + cut
        deltay = (upper-lower)/partition
        partition_coordinates = np.array([[box[0],box[1],lower,lower,box[4],box[5]]])
        delta.append(deltay)
        for i in range(partition):
            new_row = np.array([box[0],box[1],box[2],lower+(i+1)*deltay,box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append((i+1)*deltay)
        return partition_coordinates,delta
    elif direction == 'z':
        cut = cutoff*(box[5]-box[4])
        delta = []
        upper = box[5] - cut
        lower = box[4] + cut
        deltaz = (upper-lower)/partition
        partition_coordinates = np.array([[box[0],box[1],box[2],box[3],lower,lower]])
        delta.append(deltaz)
        for i in range(partition):
            new_row = np.array([box[0],box[1],box[2],box[3],lower,lower+(i+1)*deltaz])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append((i+1)*deltaz)
        return partition_coordinates, delta
    else:
        print "getpartioningbtt: Partition parameters not valid."
        


# partitioning: CENTER TO OUTSIDE (cto)
# input:
# - sample_box = (xmin, xmax, ymin, ymax, zmin, zmax)
# - direction = 'x'/'y'/'z'
# - partition = int


def getpartitioningcto(box, direction, partition, cutoff=0):
    if direction == 'x':
        cut = cutoff*(box[1]-box[0])
        delta = []
        upper = box[1] - cut
        lower = box[0] + cut
        deltax = (upper - lower)/partition/2.
        centerx = (upper + lower)/2 
        partition_coordinates = np.array([[centerx,centerx,box[2],box[3],box[4],box[5]]])
        delta.append(0)
        for i in range(partition):
            new_row = np.array([centerx-(i+1)*deltax,centerx+(i+1)*deltax,box[2],box[3],box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(2*(i+1)*deltax)
        return partition_coordinates, delta
    elif direction == 'y':
        cut = cutoff*(box[3]-box[2])
        delta = []
        upper = box[3] - cut
        lower = box[2] - cut
        deltay = (upper - lower)/partition/2.
        centery = (upper + lower)/2
        partition_coordinates = np.array([[box[0],box[1],centery,centery,box[4],box[5]]])
        delta.append(0)
        for i in range(partition):
            new_row = np.array([box[0],box[1],centery-(i+1)*deltay,centery+(i+1)*deltay,box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(2*(i+1)*deltay)
        return partition_coordinates,delta
    elif direction == 'z':
        cut = cutoff*(box[5]-box[4])
        delta = []
        upper = box[5] - cut
        lower = box[4] + cut
        deltaz = (upper - lower)/partition/2.
        centerz = (upper + lower)/2
        partition_coordinates = np.array([[box[0],box[1],box[2],box[3],centerz,centerz]])
        delta.append(0)
        for i in range(partition):
            new_row = np.array([box[0],box[1],box[2],box[3],centerz-(i+1)*deltaz,centerz+(i+1)*deltaz])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(2*(i+1)*deltaz)
        return partition_coordinates, delta
    else:
        print "getpartioningcto: Partition parameters not valid."



# partitioning: EQUIDISTAND SLICING
# input:
# - sample_box (vtk) = (xmin, xmax, ymin, ymax, zmin, zmax)
# - direction = 'x'/'y'/'z'
# - partition = int


def getpartitioningequal(box, direction, partition, cutoff=0):
    if direction == 'x':
        delta = []
        deltax = (box[1]-box[0])/partition
        lower = box[0]
        upper = box[0]+deltax
        partition_coordinates = np.array([[lower,upper,box[2],box[3],box[4],box[5]]])
        delta.append(deltax)
        for i in range(partition-1):
            new_row = np.array([lower+(i+1)*deltax,upper+(i+1)*deltax,box[2],box[3],box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(deltax)
        return partition_coordinates, delta
    elif direction == 'y':
        delta =[]
        deltay = (box[3]-box[2])/partition
        lower = box[2]
        upper = box[2]+deltay
        partition_coordinates = np.array([[box[0],box[1],lower,upper,box[4],box[5]]])
        delta.append(deltay)
        for i in range(partition-1):
            new_row = np.array([box[0],box[1],lower+(i+1)*deltay,upper+(i+1)*deltay,box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(deltay)
        return partition_coordinates,delta
    elif direction == 'z':
        delta = []
        deltaz = (box[5]-box[4])/partition
        lower = box[4]
        upper = box[4]+deltaz
        partition_coordinates = np.array([[box[0],box[1],box[2],box[3],lower,upper]])
        delta.append(deltaz)
        for i in range(partition-1):
            new_row = np.array([box[0],box[1],box[2],box[3],lower+(i+1)*deltaz,upper+(i+1)*deltaz])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(deltaz)
        return partition_coordinates, delta
    else:
        print "getpartioningequal: Partition parameters not valid."    

    

# partitioning: EQUIDISTAND SLICING (cto)
# input:
# - sample_box (vtk) = (xmin, xmax, ymin, ymax, zmin, zmax)
# - direction = 'x'/'y'/'z'
# - partition = int
#

def getpartitioningttb(box, direction, partition, cutoff=0):
    """
    Parameters:
    -----------
    box : [[xmin, xmax, ymin, ymax, zmin, zmax]]
    
    direction : string {x, y, z}
    
    partition : int
        number of desired partitions
        
    cutoff : optional float
        cuts off boundary to avoid boundary effects
    
    Return:
    -------
    partition_coordinates : np.array
        rows of: xmin, xmax, ymin, ymax, zmin, zmax
    
    delta : float
        spacing of partitioning
    """
    if direction == 'x':
        cut = cutoff*(box[1]-box[0])
        delta = []
        # cut boundaries off
        upper = round(box[1] - cut)
        lower = round(box[0] + cut)
        deltax = (upper-lower)/partition
        lower = upper-deltax
        partition_coordinates = np.array([[lower,upper,box[2],box[3],box[4],box[5]]])
        delta.append(deltax)
        for i in range(partition-1):
            new_row = np.array([lower-(i+1)*deltax,upper,box[2],box[3],box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(deltax)
        return partition_coordinates, delta
    elif direction == 'y':
        cut = cutoff*(box[3]-box[2])
        delta =[]
        upper = box[3] - cut
        lower = box[2] + cut
        deltay = (upper-lower)/partition
        lower = upper - deltay
        partition_coordinates = np.array([[box[0],box[1],lower,upper,box[4],box[5]]])
        delta.append(deltay)
        for i in range(partition-1):
            new_row = np.array([box[0],box[1],lower-(i+1)*deltay,upper,box[4],box[5]])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(deltay)
        return partition_coordinates,delta
    elif direction == 'z':
        cut = cutoff*(box[5]-box[4])
        delta = []
        upper = box[5] - cut
        lower = box[4] + cut
        deltaz = (upper-lower)/partition
        lower = upper-deltaz
        partition_coordinates = np.array([[box[0],box[1],box[2],box[3],lower,upper]])
        delta.append(deltaz)
        for i in range(partition-1):
            new_row = np.array([box[0],box[1],box[2],box[3],lower-(i+1)*deltaz,upper])
            partition_coordinates = np.vstack((partition_coordinates,new_row))
            delta.append(deltaz)
        return partition_coordinates, delta
    else:
        print "getpartioningttb: Partition parameters not valid."    
