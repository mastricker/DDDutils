# Helper functions for reading in variables from file

import fileinput

def get_dimensions(pathtoinfos):

    """
    Helper function to read out dimensions of a config from the file
    'start_infos.dat'

    pathtoinfos: absolute path to start_infos.dat file
    """

    for line in fileinput.input(pathtoinfos):

        if line.startswith('XLEN'):
            xlen = line.strip().split()[2]
            xlen = float(xlen[0:xlen.rfind(',')])
            
        elif line.startswith('YLEN'):
            ylen = line.strip().split()[2]
            ylen = float(ylen[0:ylen.rfind(',')])
            
        elif line.startswith('ZLEN'):
            zlen = line.strip().split()[2]
            zlen = float(zlen[0:zlen.rfind(',')])
    
    return xlen, ylen, zlen

