import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# read provided netcdf
ncpath='../ncfiles/Ocean1_input_geom_v1.01.nc'
x = netCDF4.Dataset(ncpath).variables['x'][::2]
y = netCDF4.Dataset(ncpath).variables['y'][::2]
[X,Y]=np.meshgrid(x/1.0e3,y/1.0e3)

# read 3D file
file3D = netCDF4.Dataset('../ncfiles/Ideal_3D.nc','r+')
area3D=file3D.variables['area'][:]
# read 3D file
file2D = netCDF4.Dataset('../ncfiles/Ideal_2D.nc','r+')
area2D=file2D.variables['area'][:]

# read ocean geometry
depth = netCDF4.Dataset('../ncfiles/ocean_geometry.nc').variables['D'][:]

xx=np.array([430,450,700])*1.0e3
yy=[-760,-750,-40]

p = np.poly1d(np.polyfit(xx, yy, 1))
lower=p(x)
plt.figure()
plt.plot(x/1.0e3,lower)
plt.plot(xx/1.0e3,yy,'ro')
plt.xlim(320,800)
plt.show()
upper=np.ones(len(x))*10.
a=np.nonzero(x<=430*1.0e3)[-1][-1]
lower[0:a+1]=-760;
a=np.nonzero(x>=700*1.0e3)[-1][0]
lower[a::]=-40
 
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
plt.plot(x/1.0e3,-depth[:,0])
#plt.plot(x/1.0e3,lower3D[:,0],'r')
#plt.plot(x/1.0e3,lower3D[:,20],'b')
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

