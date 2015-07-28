"""
Plotting library for DDD raw data.
"""

# relative imports
from .. import constants
from .  import util_plot as up


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc, pyplot, cm
import numpy as np
import os,sys
from subprocess import Popen, PIPE

import util_plot as up

# Utils for glide system plotting fcc

gsystems = ("A2", "A3", "A6",\
            "B2", "B4", "B5",\
            "C1", "C3", "C5",\
            "D1", "D4", "D6")

map_schmidboas = (2,3,1,6,4,5,7,9,8,11,10,12)

##### plot
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=nsystem)
cNorm = matplotlib.colors.Normalize(0,nsystem)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

plt.rc('lines', linewidth=1)
rc('text', usetex=True)
rc('text.latex', preamble = ','.join('''
\usepackage{txfonts}
\usepackage{lmodern}
'''.split()))
rc('font', family='sans-serif', weight='normal', style='normal')
