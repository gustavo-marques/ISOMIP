# All the plotting related functions used in the notebooks are kept here
# author: Gustavo M. Marques

from matplotlib import pyplot as plt
import sys
sys.path.append('~/python/pyGVtools/')
import m6toolbox
import netCDF4
import subprocess
import numpy as np

# Define a function to plot a section
def plot_section(file_handle, record, xq, i=0, variable='salt',eta='e',clim=(33.8,34.55), plot_grid=True, rep='pcm', xlim=(320,800), ylim=(-720,0), cmap=plt.cm.jet):
    """Plots a section of by reading vertical grid and scalar variable and super-sampling
    both in order to plot vertical and horizontal reconstructions.
    
    Optional arguments have defaults for plotting salinity and overlaying the grid.
    """
    e = file_handle.variables[eta][record,:,:,i] # Vertical grid positions
    s = file_handle.variables[variable][record,:,:,i] # Scalar field to color
    x,z,q = m6toolbox.section2quadmesh(xq, e, s, representation=rep) # This yields three areas at twice the model resolution
    plt.pcolormesh(x, z, q, cmap=cmap);
    plt.colorbar()
    plt.clim(clim)
    if plot_grid: plt.plot(x, z.T, 'k', hold=True);
    plt.ylim(ylim)
    plt.xlim(xlim)

def plot_diags(path,n1,n2,n3,n4):
    """
    Plot some diagnostics computed by the model. These are: total mass, energy,
     ave. temp/salt, and maxCFL as a function of time.
    """
    eval(str("subprocess.call(['./get_output.csh %s'])"% (path+'layer/'+n1)))

    subprocess.call(['./get_output.csh' + path])

def plot_stats(file_handle,variable='v',labels=['layer','rho','sigma','z']):
    # Two subplots, the axes array is 1-d
    plt.figure(figsize=(18,9))
    plt.title('Max. absolute vel.')
    units = file_handle[0].variables['u'].units
    for i in range(len(labels)):
        time = file_handle[i].variables['Time'][:]
        smin=np.zeros(len(time)); smax=np.zeros(len(time))
        total=np.zeros(len(time))
        s = file_handle[i].variables[variable][:]      
        for t in range(len(time)):
            smin[t]=s[t,:].min(); smax[t]=s[t,:].max()
            if smax[t]>=np.abs(smin[t]):
               total[t]=smax[t]
               plt.semilogy(time[t],total[t],'+',ms=2)
            else:
               total[t]=np.abs(smin[t])
               plt.semilogy(time[t],total[t],'o',ms=2)

        plt.plot(time,total,lw=2.0, label=labels[i])
        plt.legend(loc=0, shadow=True)
        plt.xlabel('Time') 
        plt.ylabel(units) 

def get_rho(S,T,rho0=1027.51,Sref=34.2,Tref=-1.0,alpha = 3.733e-5, beta = 7.843e-4):
    """ Compute linear EoS """
    rho = rho0 * (1 - alpha * (T-Tref) + beta * (S-Sref))
    return rho
