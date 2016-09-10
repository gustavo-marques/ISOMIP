import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

# read provided netcdf
ncpath='../ncfiles/Ocean1_input_geom_v1.01.nc'
x = netCDF4.Dataset(ncpath).variables['x'][::2]
y = netCDF4.Dataset(ncpath).variables['y'][::2]
[X,Y]=np.meshgrid(x/1.0e3,y/1.0e3)

# read 3D file
file3D = netCDF4.Dataset('../ncfiles/Ideal_step_3D.nc','r+')
area3D=file3D.variables['area'][:]
# read 3D file
file2D = netCDF4.Dataset('../ncfiles/Ideal_step_2D.nc','r+')
area2D=file2D.variables['area'][:]


upper=np.zeros(len(x));lower=np.zeros(len(x))
a=np.nonzero(x<=700*1.0e3)[-1][-1]
Dmax=-300
lower[0:a]=Dmax;
#lower[0:a-5]=Dmax;
#lower[a-5]=Dmax*0.875;lower[a-4]=Dmax*3/4.;lower[a-3]=Dmax/2.;lower[a-2]=Dmax/4.;lower[a-1]=Dmax/8.;lower[a]=-2.0;
upper[0:a]=10;

#apply filter
sigma = [6] #  the standard deviation of the distribution
lower = gaussian_filter(lower,sigma)
upper = gaussian_filter(upper,sigma)
 
thick = upper - lower
# plot
plt.figure(figsize=(12,6))
plt.plot(x/1.0e3,lower,'g')
plt.show()

# put into arays
area2D = np.ones((area2D.shape))* (x[1]-x[0]) * (y[1]-y[0])
area3D = np.ones((area3D.shape))* (x[1]-x[0]) * (y[1]-y[0])
upper2D = np.ones((area2D.shape))
upper3D = np.ones((area3D.shape))
lower2D = np.ones((area2D.shape))
lower3D = np.ones((area3D.shape))

# 3D
for i in range(area3D.shape[1]):
        upper3D[:,i]=upper[:]
        lower3D[:,i]=lower[:]
    
thick3D=upper3D - lower3D
area3D[thick3D==0.0]=0.0

# 2D
for i in range(area2D.shape[1]):
    upper2D[:,i]=upper[:]
    lower2D[:,i]=lower[:]
    
thick2D=upper2D - lower2D
area2D[thick2D==0.0]=0.0

plt.figure()
plt.plot(x/1.0e3,lower2D[:,0])
plt.show()

# put into netcdf files

# 3D
file3D.variables['area'][:,:] = area3D[:,:]
file3D.variables['thick'][:,:] = thick3D[:,:]


# 2D (middle of the domain)
file2D.variables['area'][:,:] = area2D[:,:]
file2D.variables['thick'][:,:] = thick2D[:,:]


file2D.close()
file3D.close()

