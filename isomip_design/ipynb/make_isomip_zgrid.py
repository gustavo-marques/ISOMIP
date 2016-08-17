# Setup the vertical grid to be used in the interpolation of the ISOMIP+ simulations 
# Gustavo Marques, Aug 2016

import netCDF4
import numpy as np

# read ncfile
file = netCDF4.Dataset('../ncfiles/ISOMIP_zgrid.nc','r+')
file.variables['zt'][:] = np.arange(2.5,720.,5.)
file.variables['zw'][:] = np.arange(0,725.,5.)
file.close()
