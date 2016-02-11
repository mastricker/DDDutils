# Helper functions for reading in variables from file

import fileinput

def get_dimensions(pathtoinfos):

    """
    Helper function to read out dimensions of a config from the file
    'start_infos.dat'

    pathtoinfos: absolute path to start_infos.dat file
    """

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
