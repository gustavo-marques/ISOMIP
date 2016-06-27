# All the plotting related functions used in the notebooks are kept here
# author: Gustavo M. Marques

from matplotlib import pyplot as plt
import sys
sys.path.append('~/python/pyGVtools/')
import m6toolbox
import netCDF4

# Define a function to plot a section
def plot_section(file_handle, record, xq, variable='salt', clim=(33.8,34.55), plot_grid=True, rep='linear', xlim=(400,800), ylim=(-720,0), cmap=plt.cm.jet):
    """Plots a section of by reading vertical grid and scalar variable and super-sampling
    both in order to plot vertical and horizontal reconstructions.
    
    Optional arguments have defaults for plotting salinity and overlaying the grid.
    """
    e = file_handle.variables['e'][record,:,0,:] # Vertical grid positions
    s = file_handle.variables[variable][record,:,0,:] # Scalar field to color
    x,z,q = m6toolbox.section2quadmesh(xq, e, s, representation=rep) # This yields three areas at twice the model resolution
    plt.pcolormesh(x, z, q, cmap=cmap);
    plt.colorbar()
    plt.clim(clim)
    if plot_grid: plt.plot(x, z.T, 'k', hold=True);
    plt.ylim(ylim)
    plt.xlim(xlim)



