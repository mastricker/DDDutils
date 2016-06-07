# Helper functions for reading in variables from file

import fileinput

def get_dimensions(pathtoinfos):

    """
    Helper function to read out dimensions of a config from the file
    'start_infos.dat'

    pathtoinfos: absolute path to start_infos.dat file
    """

    print '---', pathtoinfos[-15:]

    if pathtoinfos[-15:] != 'start_infos.dat':
        if pathtoinfos[-16] == '/':
            pathtoinfos += 'start_infos.dat'
        else:
            pathtoinfos += '/start_infos.dat'
        
    for line in fileinput.input(pathtoinfos):

        if line.lstrip().startswith('XLEN'):
            xlen = line.strip().split()[2]
            xlen = float(xlen[0:xlen.rfind(',')])
            
        elif line.lstrip().startswith('YLEN'):
            ylen = line.strip().split()[2]
            ylen = float(ylen[0:ylen.rfind(',')])
            
        elif line.lstrip().startswith('ZLEN'):
            zlen = line.strip().split()[2]
            zlen = float(zlen[0:zlen.rfind(',')])

    return xlen, ylen, zlen

def get_material_parameters(pathtoinfos):

    """
    Helper function to read out material parameter of a config from the file
    'start_infos.dat'

    pathtoinfos: absolute path to start_infos.dat file
    """

    pathtoinfos = pathtoinfos + '/start_infos.dat'

    print 'Reading material param from', pathtoinfos

    for line in fileinput.input(pathtoinfos):

        if line.lstrip().startswith('MU'):
            mu = float( line.strip().split()[2].replace(',','') )

        elif line.lstrip().startswith('NU'):
            nu = float( line.strip().split()[2].replace(',','') )
 
        elif line.lstrip().startswith('LATTICE_CONSTANT'):
            a = float( line.strip().split()[2].replace(',','') )

    return mu, nu, a

def get_ids_from_simulations(pathtocfg):
    """
    Helper function, which returns the correct indices for a given path.
    
    pathtocfg : absolute path to simulation config folder
    """

    dict_init = {'directstrain' : 0, 'frsources' : 1, 'restart7500' : 2}
    dict_size = {'AR3_1' : 1, 'AR3_2' : 2, 'AR3_3' : 3, 'AR3_4' : 4}


    p = pathtocfg.split('/')

    init = p[5]
    size = p[6]
    

    id_init = dict_init[init]
    id_size = dict_size[size]
    id_dir = int(p[7])
    
    return id_init, id_size, id_dir
