import netCDF4
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt

# read provided netcdf
ncpath='../ncfiles/Ocean1_input_geom_v1.01.nc'
x = netCDF4.Dataset(ncpath).variables['x'][:]
y = netCDF4.Dataset(ncpath).variables['y'][:]
upperSurface = netCDF4.Dataset(ncpath).variables['upperSurface'][:]
lowerSurface = netCDF4.Dataset(ncpath).variables['lowerSurface'][:]
[X,Y]=np.meshgrid(x/1.0e3,y/1.0e3)

# read 3D file
file3D = netCDF4.Dataset('../ncfiles/Ocean1_3D.nc','r+')
# read 3D file
file2D = netCDF4.Dataset('../ncfiles/Ocean1_2D.nc','r+')

thick = upperSurface - lowerSurface

lowerSurface[thick<100.]=0.0
upperSurface[thick<100.]=0.0
thick = upperSurface - lowerSurface

# original data
plt.figure(figsize=(12,6))
plt.contourf(X,Y,upperSurface,50)
plt.colorbar()
plt.title('upperSurface [m]')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,thick,50)
plt.colorbar()
plt.title('Thickness [m]')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,lowerSurface, 50)
plt.colorbar()
plt.title('lowerSurface [m]')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
#plt.show()

# find ice shelf front, it does not depend on y
front=np.nonzero(lowerSurface[0,:]==0)[0][0]
#front=np.nonzero(lowerSurface[0,:]==0)[0][-1]

# smoothing
lowerSurface_smoth = np.zeros(lowerSurface.shape)
upperSurface_smoth = np.zeros(upperSurface.shape)

sigma = [4,4] #  the standard deviation of the distribution
#lowerSurface_smoth[:,0:front] = gaussian_filter(lowerSurface[:,0:front],sigma)
lowerSurface_smoth = gaussian_filter(lowerSurface,sigma)
#upperSurface_smoth[:,0:front] = gaussian_filter(upperSurface[:,0:front],sigma)
upperSurface_smoth = gaussian_filter(upperSurface,sigma)

thick_smoth = upperSurface_smoth - lowerSurface_smoth

# plot after first smoothing
plt.figure(figsize=(12,6))
plt.contourf(X,Y,upperSurface_smoth - upperSurface,50)
plt.colorbar()
plt.title('upperSurface (differece)')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,thick_smoth - thick,50)
plt.colorbar()
plt.title('Thickness (differece)')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,lowerSurface_smoth - lowerSurface, 50)
plt.colorbar()
plt.title('lowerSurface (difference)')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
#plt.show()

# smooth the shelf fron, just for layer model
#sigma = 1 # this seems to be reasonable
#front10=front-10
#for j in range(lowerSurface.shape[0]):
#    lowerSurface_smoth[j,front+1]=-50;lowerSurface_smoth[j,front+2]=-25.;lowerSurface_smoth[j,front+3]=-12.5
    #lowerSurface_smoth[j,front10::] = gaussian_filter(lowerSurface_smoth[j,front10::],sigma)
    #upperSurface_smoth[j,front::] = gaussian_filter(upperSurface_smoth[j,front::],sigma)

    
#thick_smoth = upperSurface_smoth - lowerSurface_smoth
#lowerSurface_smoth[thick<100.]=0.0
#thick_smoth = upperSurface_smoth - lowerSurface_smoth

# plot again
plt.figure(figsize=(12,6))
plt.contourf(X,Y,upperSurface_smoth - upperSurface,50)
plt.colorbar()
plt.title('upperSurface (differece)')
plt.xlabel('x [km]')
plt.ylabel('y [km]')

plt.figure(figsize=(12,6))
plt.contourf(X,Y,lowerSurface_smoth - lowerSurface, 50)
plt.colorbar()
plt.title('lowerSurface (difference)')
plt.xlabel('x [km]')
plt.ylabel('y [km]')

plt.figure(figsize=(12,6))
plt.contourf(X,Y,thick_smoth - thick,50)
plt.colorbar()
plt.title('Thickness (differece)')
plt.close('all')

plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[40,:],'b');plt.plot(lowerSurface[40,:],'r')
plt.xlim(300,350)
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at x = 40 km')

plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[5,:],'b');plt.plot(lowerSurface[5,:],'r')
plt.xlim(300,350)
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at x = 5 km')

plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[:,100],'b');plt.plot(lowerSurface[:,100],'r')
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at y = 400 km')


plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[:,200],'b');plt.plot(lowerSurface[:,200],'r')
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at y = 500 km')
plt.show()

# put into coarse grid
x1=x[::2];y1=y[::2]
thick1=thick_smoth[::2,::2]
height=lowerSurface_smoth[::2,::2]
top=upperSurface_smoth[::2,::2]
height[thick1<100.]=0.0
top[thick1<100.]=0.0
thick1 = top - height

sigma = [2,2] #  the standard deviation of the distribution
thick1 = gaussian_filter(thick1,sigma)
height = gaussian_filter(height,sigma)

height[thick1<15.]=0.0
top[thick1<15.]=0.0
thick1 = top - height

area = np.ones((thick1.shape))* (x1[1]-x1[0]) * (y1[1]-y1[0])
area[thick1==0.0]=0.0

#save into netcdf file
# 3D
file3D.variables['area'][:,:] = area[:,:].T
file3D.variables['thick'][:,:] = thick1[:,:].T
file3D.variables['height'][:,:] = height[:,:].T

# 2D (middle of the domain)
for i in range(file2D.variables['area'].shape[1]):
    file2D.variables['area'][:,i] = area[20,:]
    file2D.variables['thick'][:,i] = thick1[20,:]
    file2D.variables['height'][:,i] = height[20,:] 

file2D.close()
file3D.close()
