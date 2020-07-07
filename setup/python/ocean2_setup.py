import netCDF4
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from scipy import interpolate


# ice shelf
# read provided netcdf
ncpath='../ncfiles/Ocean2_input_geom_v1.01.nc'

# read IC files
file3Dic = netCDF4.Dataset('../ncfiles/Ocean2_3D.nc','r+')

x = netCDF4.Dataset(ncpath).variables['x'][:]
y = netCDF4.Dataset(ncpath).variables['y'][:]
[X,Y]=np.meshgrid(x/1.0e3,y/1.0e3)

upperSurface = netCDF4.Dataset(ncpath).variables['upperSurface'][:]
lowerSurface = netCDF4.Dataset(ncpath).variables['lowerSurface'][:]

thick = upperSurface - lowerSurface

# smoothing
lowerSurface_smoth = np.zeros(lowerSurface.shape)
upperSurface_smoth = np.zeros(upperSurface.shape)

sigma = [2,2] # (y,x) the standard deviation of the distribution
#lowerSurface_smoth[:,0:front] = gaussian_filter(lowerSurface[:,0:front],sigma)
lowerSurface_smoth = gaussian_filter(lowerSurface,sigma)
#upperSurface_smoth[:,0:front] = gaussian_filter(upperSurface[:,0:front],sigma)
upperSurface_smoth = gaussian_filter(upperSurface,sigma)

thick_smoth = upperSurface_smoth - lowerSurface_smoth

# calve
thick_smoth[thick_smoth<100] = 0.0

print(thick_smoth.shape)

# interpolate to coarse grid
xnew=x[::2];ynew=y[::2]
f_thick = interpolate.interp2d(x, y, thick_smoth, kind='cubic')
thick_new=f_thick(xnew, ynew)

# simple calving
thick_new[thick_new<100.] = 0.0
thick_new = gaussian_filter(thick_new,sigma)
# calving again to remove thin ice shelf
thick_new[thick_new<100.] = 0.0

x2=xnew/1.0e3;y2=ynew/1.0e3
jm,im = thick_new.shape

# compute mass (kg/m^2)
rho_ice = 918.
mass = thick_new * rho_ice
# update area
area = np.ones((thick_new.shape))* (xnew[1]-xnew[0]) * (ynew[1]-ynew[0])
area[thick_new==0]=0.0

#save into netcdf files
print('Saving initial conditions...')
file3Dic.variables['area'][:] = area[:,:]
file3Dic.variables['thick'][:] = thick_new[:,:]

print('Done!')

file3Dic.close()
